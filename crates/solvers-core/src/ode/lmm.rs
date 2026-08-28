//! The linear multistep stepper.
//!
//! The coefficients are not stored, they are solved for at every step from the
//! order conditions on the actual step size history. That makes the variable
//! step behaviour a property of the framework rather than of the individual
//! method, and it is what lets BDF, Adams-Bashforth, Adams-Moulton and any
//! other pattern of free coefficients share this one implementation.
//!
//! The order ramps up as the history fills, so no separate start up method is
//! needed: a k-step family run with one stored point is its own one-step
//! member.

use super::newton_matrix::NewtonMatrix;
use super::rk::RkStepper;
use super::{Options, StepOutcome, Stepper};
use crate::method::{LmmCoefficients, LmmFamily};
use crate::nonlinear::{NonlinearSolver, Residual, SolverKind};
use crate::problem::{Problem, Stats};
use crate::simd;
use std::collections::VecDeque;

pub struct LmmStepper {
    family: LmmFamily,
    dim: usize,
    rtol: f64,
    atol: f64,
    max_jacobian_age: u32,

    y: Vec<f64>,
    y_new: Vec<f64>,
    f_new: Vec<f64>,
    /// Past solution values, most recent first.
    y_history: VecDeque<Vec<f64>>,
    /// Past derivative values, most recent first.
    f_history: VecDeque<Vec<f64>>,
    /// Past step sizes, most recent first.
    h_history: VecDeque<f64>,

    scale: Vec<f64>,
    rhs_const: Vec<f64>,
    work: Vec<f64>,
    low: Vec<f64>,

    linear: NewtonMatrix,
    nonlinear: NonlinearSolver,

    /// One step method used to fill the history for families that cannot ramp
    /// their own order down, such as Nystrom and Milne-Simpson.
    startup: Option<RkStepper>,
    /// True while the last attempted step came from the start up method.
    in_startup: bool,

    /// Number of steps the formula currently uses.
    order_in_use: usize,
    /// Order suggested for the next step by the error comparison.
    next_order: usize,
    step_sizes: Vec<f64>,
    last_h: f64,
    have_step: bool,
}

impl LmmStepper {
    pub fn new(family: &LmmFamily, dim: usize, options: &Options) -> LmmStepper {
        LmmStepper::with_startup(family, dim, options, None)
    }

    pub fn with_startup(
        family: &LmmFamily,
        dim: usize,
        options: &Options,
        startup: Option<RkStepper>,
    ) -> LmmStepper {
        let k = family.steps;
        LmmStepper {
            family: family.clone(),
            dim,
            rtol: options.rtol,
            atol: options.atol,
            max_jacobian_age: options.max_jacobian_age,
            y: vec![0.0; dim],
            y_new: vec![0.0; dim],
            f_new: vec![0.0; dim],
            y_history: VecDeque::with_capacity(k + 1),
            f_history: VecDeque::with_capacity(k + 1),
            h_history: VecDeque::with_capacity(k + 1),
            scale: vec![1.0; dim],
            rhs_const: vec![0.0; dim],
            work: vec![0.0; dim],
            low: vec![0.0; dim],
            linear: NewtonMatrix::new(dim),
            nonlinear: NonlinearSolver::new(options.nonlinear),
            startup,
            in_startup: false,
            order_in_use: 1,
            next_order: 1,
            step_sizes: vec![0.0; k + 1],
            last_h: 0.0,
            have_step: false,
        }
    }

    /// Steps of history that are actually available.
    fn available(&self) -> usize {
        self.y_history.len()
    }

    fn fill_step_sizes(&mut self, h: f64, k: usize) {
        self.step_sizes.resize(k + 1, h.abs());
        self.step_sizes[0] = h.abs();
        for j in 1..=k {
            let index = j - 1;
            self.step_sizes[j] = if index < self.h_history.len() {
                self.h_history[index].abs()
            } else {
                h.abs()
            };
        }
    }

    /// `rhs_const = h * sum_{j>=1} beta_j f_{n-j} - sum_{j>=1} alpha_j y_{n-j}`
    fn assemble(&mut self, coefficients: &LmmCoefficients, h: f64, out: &mut [f64]) {
        let k = coefficients.alpha.len() - 1;
        for value in out.iter_mut() {
            *value = 0.0;
        }
        for j in 1..=k {
            let index = j - 1;
            if index >= self.y_history.len() {
                break;
            }
            let alpha = coefficients.alpha[j];
            let beta = coefficients.beta[j];
            let y_past = &self.y_history[index];
            let f_past = &self.f_history[index];
            for i in 0..out.len() {
                out[i] += h * beta * f_past[i] - alpha * y_past[i];
            }
        }
    }
}

impl<P: Problem + ?Sized> Stepper<P> for LmmStepper {
    fn control_order(&self) -> usize {
        self.order_in_use
    }

    fn state(&self) -> &[f64] {
        &self.y
    }

    fn proposed(&self) -> &[f64] {
        &self.y_new
    }

    fn is_adaptive(&self) -> bool {
        self.family.steps > 1
    }

    fn max_growth(&self) -> f64 {
        // The coefficients stay well conditioned only if the history spacing
        // changes gradually.
        2.0
    }

    /// Shrink the step while the formula is still ramping up its order.
    ///
    /// A k-step method has no history to start from, so the first steps run at
    /// reduced order q and would contribute a local error of order q + 1 to a
    /// solution that is meant to be order p. Taking those steps with
    /// `h^((p+1)/(q+1))` instead puts their local error at the same level as
    /// the full order method's, which preserves the global order without
    /// needing a second method to start with. The variable step coefficients
    /// then absorb the uneven history spacing on their own.
    fn step_limit(&self, h: f64) -> f64 {
        // A dedicated start up method already produces the right local error,
        // so there is nothing to compensate for.
        if self.startup.is_some() {
            return h;
        }
        let target = self.family.steps;
        let current = self.available().min(target);
        if current >= target || h <= 0.0 {
            return h;
        }
        let magnitude = h.abs();
        if magnitude >= 1.0 {
            return h;
        }
        let exponent = (target + 1) as f64 / (current + 1) as f64;
        let reduced = magnitude.powf(exponent).max(1e-13);
        h.signum() * reduced.min(magnitude)
    }

    fn start(&mut self, problem: &P, stats: &mut Stats, t: f64, y: &[f64]) {
        self.y.copy_from_slice(y);
        self.y_new.copy_from_slice(y);
        self.y_history.clear();
        self.f_history.clear();
        self.h_history.clear();
        self.linear.invalidate();
        self.nonlinear.reset();
        self.order_in_use = 1;
        self.next_order = 1;
        self.have_step = false;

        let mut f0 = vec![0.0; self.dim];
        stats.rhs_evals += 1;
        problem.rhs(t, y, &mut f0);
        self.y_history.push_front(y.to_vec());
        self.f_history.push_front(f0);

        self.in_startup = false;
        if let Some(rk) = &mut self.startup {
            <RkStepper as Stepper<P>>::start(rk, problem, stats, t, y);
        }
    }

    fn attempt(&mut self, problem: &P, stats: &mut Stats, t: f64, h: f64) -> StepOutcome {
        // Fill the history with the dedicated start up method when the family
        // has no lower order member to fall back on.
        if self.startup.is_some() && self.available() < self.family.min_steps {
            self.in_startup = true;
            self.order_in_use = self.family.min_steps;
            let rk = self.startup.as_mut().expect("start up method present");
            let outcome = <RkStepper as Stepper<P>>::attempt(rk, problem, stats, t, h);
            if !outcome.ok {
                return outcome;
            }
            self.y_new
                .copy_from_slice(<RkStepper as Stepper<P>>::proposed(
                    self.startup.as_ref().expect("start up method present"),
                ));
            stats.rhs_evals += 1;
            let mut f_new = std::mem::take(&mut self.f_new);
            problem.rhs(t + h, &self.y_new, &mut f_new);
            self.f_new = f_new;
            self.last_h = h;
            return outcome;
        }
        self.in_startup = false;

        let k = self
            .next_order
            .min(self.family.steps)
            .min(self.available())
            .max(1);
        self.order_in_use = k;

        let Some(family) = self.family.with_steps(k) else {
            return StepOutcome::failed();
        };
        self.fill_step_sizes(h, k);
        let Ok(coefficients) = family.coefficients(&self.step_sizes) else {
            return StepOutcome::failed();
        };

        let mut rhs_const = std::mem::take(&mut self.rhs_const);
        self.assemble(&coefficients, h, &mut rhs_const);

        let t_new = t + h;
        let beta0 = coefficients.beta[0];

        if beta0.abs() < 1e-15 {
            // Explicit formula, the update is already assembled.
            self.y_new.copy_from_slice(&rhs_const);
            stats.rhs_evals += 1;
            let mut f_new = std::mem::take(&mut self.f_new);
            problem.rhs(t_new, &self.y_new, &mut f_new);
            self.f_new = f_new;
        } else {
            // Predict with an explicit Euler step off the most recent point.
            let last_y = &self.y_history[0];
            let last_f = &self.f_history[0];
            for i in 0..self.dim {
                self.y_new[i] = last_y[i] + h * last_f[i];
            }

            simd::error_scale(self.atol, self.rtol, &self.y, &self.y, &mut self.scale);
            let full_newton = self.nonlinear.config.kind == SolverKind::Newton;

            let mut residual = LmmResidual {
                problem,
                stats,
                linear: &mut self.linear,
                t: t_new,
                h_beta0: h * beta0,
                rhs_const,
                work: std::mem::take(&mut self.work),
                max_age: self.max_jacobian_age,
                full_newton,
            };

            let mut candidate = std::mem::take(&mut self.y_new);
            self.nonlinear.reset();
            let outcome = self.nonlinear.solve(&mut residual, &mut candidate, &self.scale);

            rhs_const = residual.rhs_const;
            self.work = residual.work;
            self.y_new = candidate;
            outcome.record(stats);

            if !outcome.converged() {
                self.rhs_const = rhs_const;
                return StepOutcome::failed();
            }
            stats.rhs_evals += 1;
            let mut f_new = std::mem::take(&mut self.f_new);
            problem.rhs(t_new, &self.y_new, &mut f_new);
            self.f_new = f_new;
        }
        self.rhs_const = rhs_const;

        if !self.y_new.iter().all(|v| v.is_finite()) {
            return StepOutcome::failed();
        }

        // Error estimate: apply the order reduced formula to the converged
        // derivative. That costs no extra solve and has the right asymptotics.
        let mut error = 0.0;
        if k > 1 {
            if let Some(lower) = self.family.with_steps(k - 1) {
                self.fill_step_sizes(h, k - 1);
                if let Ok(low_coefficients) = lower.coefficients(&self.step_sizes) {
                    let mut low = std::mem::take(&mut self.low);
                    self.assemble(&low_coefficients, h, &mut low);
                    for i in 0..self.dim {
                        low[i] += h * low_coefficients.beta[0] * self.f_new[i];
                        low[i] = self.y_new[i] - low[i];
                    }
                    simd::error_scale(self.atol, self.rtol, &self.y, &self.y_new, &mut self.scale);
                    error = simd::weighted_rms(&low, &self.scale);
                    self.low = low;
                }
            }
        }

        self.last_h = h;
        StepOutcome { ok: true, error }
    }

    fn commit(&mut self, t: f64, h: f64) {
        if self.in_startup {
            // The starter is only alive until the history is full, so it needs
            // to follow the solution during that phase and not after.
            if let Some(rk) = &mut self.startup {
                <RkStepper as Stepper<P>>::commit(rk, t, h);
            }
        }
        self.y.copy_from_slice(&self.y_new);
        self.y_history.push_front(self.y_new.clone());
        self.f_history.push_front(self.f_new.clone());
        self.h_history.push_front(h);
        while self.y_history.len() > self.family.steps + 1 {
            self.y_history.pop_back();
            self.f_history.pop_back();
        }
        while self.h_history.len() > self.family.steps + 1 {
            self.h_history.pop_back();
        }
        self.linear.advance_age();
        self.have_step = true;
        // Raise the order as soon as the history supports it.
        self.next_order = (self.order_in_use + 1).min(self.family.steps);
    }

    fn reject(&mut self, h: f64) {
        if self.in_startup {
            if let Some(rk) = &mut self.startup {
                <RkStepper as Stepper<P>>::reject(rk, h);
            }
            return;
        }
        self.nonlinear.reset();
        // A rejected step at high order is often better served one order down.
        if self.order_in_use > 1 {
            self.next_order = self.order_in_use - 1;
        }
    }

    fn interpolate(&self, theta: f64, out: &mut [f64]) -> bool {
        if !self.have_step || self.f_history.is_empty() {
            return false;
        }
        // Cubic Hermite through the endpoints of the last step.
        let y0 = &self.y_history[0];
        let f0 = &self.f_history[0];
        let t2 = theta * theta;
        let t3 = t2 * theta;
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + theta;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;
        for i in 0..out.len() {
            out[i] = h00 * y0[i]
                + h10 * self.last_h * f0[i]
                + h01 * self.y_new[i]
                + h11 * self.last_h * self.f_new[i];
        }
        true
    }
}

/// `G(y) = y - h*beta_0*f(t, y) - rhs_const`
struct LmmResidual<'a, P: Problem + ?Sized> {
    problem: &'a P,
    stats: &'a mut Stats,
    linear: &'a mut NewtonMatrix,
    t: f64,
    h_beta0: f64,
    rhs_const: Vec<f64>,
    work: Vec<f64>,
    max_age: u32,
    full_newton: bool,
}

impl<'a, P: Problem + ?Sized> Residual for LmmResidual<'a, P> {
    fn dim(&self) -> usize {
        self.rhs_const.len()
    }

    fn eval(&mut self, y: &[f64], r: &mut [f64]) {
        self.stats.rhs_evals += 1;
        self.problem.rhs(self.t, y, &mut self.work);
        for i in 0..r.len() {
            r[i] = y[i] - self.h_beta0 * self.work[i] - self.rhs_const[i];
        }
    }

    fn factor(&mut self, y: &[f64]) -> bool {
        if self.full_newton {
            self.linear.invalidate();
        }
        self.linear.prepare(
            self.problem,
            self.stats,
            self.t,
            y,
            self.h_beta0,
            if self.full_newton { 0 } else { self.max_age },
        )
    }

    fn solve(&mut self, rhs: &mut [f64]) -> bool {
        self.linear.solve(rhs)
    }

    fn refresh(&mut self) {
        self.linear.invalidate();
    }
}
