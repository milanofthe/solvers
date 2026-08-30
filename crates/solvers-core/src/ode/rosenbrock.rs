//! The Rosenbrock stepper.
//!
//! What separates this from the implicit Runge-Kutta stepper next to it is that
//! there is no iteration to control. Each stage is one linear solve, decided in
//! advance, and the Jacobian is part of the formula rather than a device for
//! converging to a solution of it. A step therefore has a fixed cost: `s` right
//! hand side evaluations, one Jacobian, one factorization while the diagonal is
//! constant, and `s` back substitutions. Nothing can fail to converge, so the
//! only reason a step is rejected is that the error estimate says so.
//!
//! The form solved here is the substituted one,
//!
//! ```text
//! (I/(h gamma_i) - J) S_i = f(t + c_i h, y + sum_j a_ij S_j)
//!                           + sum_j (C_ij / h) S_j + h d_i ft
//! y1 = y + sum_i m_i S_i
//! ```
//!
//! multiplied through by `h gamma_i` so the matrix is `I - h gamma_i J`, which
//! is the same one every implicit method in this crate factors and therefore
//! shares the Jacobian bookkeeping with them.
//!
//! `ft` is `df/dt` at the start of the step, and it is not optional: without it
//! a method loses order on any problem that depends on `t` directly. It costs
//! one extra evaluation per step and only when some `d_i` is nonzero.

use super::newton_matrix::NewtonMatrix;
use super::{StepOutcome, Stepper};
use crate::method::{RosenbrockRuntime, RosenbrockTableau};
use crate::ode::Options;
use crate::problem::{Problem, Stats};
use crate::simd;

pub struct RosenbrockStepper {
    tableau: RosenbrockRuntime,
    control_order: usize,
    rtol: f64,
    atol: f64,

    y: Vec<f64>,
    y_new: Vec<f64>,
    /// The substituted stage increments `S_i`.
    s: Vec<Vec<f64>>,
    argument: Vec<f64>,
    rhs: Vec<f64>,
    /// `df/dt` at the start of the step.
    ft: Vec<f64>,
    f0: Vec<f64>,
    err: Vec<f64>,
    scale: Vec<f64>,

    linear: NewtonMatrix,
    /// Diagonal of `gamma`, one entry per stage. Every published method has a
    /// constant one and factors once, but nothing here needs that to be true.
    diagonal: Vec<f64>,
    needs_time_derivative: bool,
    last_h: f64,
}

impl RosenbrockStepper {
    pub fn new(
        tableau: &RosenbrockTableau,
        declared_order: usize,
        embedded_order: Option<usize>,
        dim: usize,
        options: &Options,
    ) -> RosenbrockStepper {
        // `options.max_jacobian_age` is deliberately not read. See `attempt`.
        let runtime = tableau.runtime();
        let stages = runtime.stages;
        // One value for the whole step where the method has one. Reading the
        // diagonal back entry by entry would hand the factorization cache `s`
        // numbers that agree to fifteen digits and differ in the sixteenth,
        // and it would then factor once per stage for no reason.
        let diagonal = match tableau.diagonal {
            Some(shared) => vec![shared.value(); stages],
            None => (0..stages).map(|i| tableau.gamma[(i, i)].value()).collect(),
        };
        let needs_time_derivative = runtime.d.iter().any(|v| v.abs() > 1e-14);

        let control_order = match embedded_order {
            Some(m) => declared_order.min(m).max(1),
            None => declared_order.max(1),
        };

        RosenbrockStepper {
            tableau: runtime,
            control_order,
            rtol: options.rtol,
            atol: options.atol,
            y: vec![0.0; dim],
            y_new: vec![0.0; dim],
            s: vec![vec![0.0; dim]; stages],
            argument: vec![0.0; dim],
            rhs: vec![0.0; dim],
            ft: vec![0.0; dim],
            f0: vec![0.0; dim],
            err: vec![0.0; dim],
            scale: vec![1.0; dim],
            linear: NewtonMatrix::new(dim),
            diagonal,
            needs_time_derivative,
            last_h: 0.0,
        }
    }

    /// `df/dt` at the start of the step, by a one sided difference in `t`.
    ///
    /// The problem interface offers a Jacobian in `y` and nothing in `t`, and a
    /// difference is what every published implementation uses here anyway. Two
    /// evaluations, once per step, not once per stage.
    fn time_derivative<P: Problem + ?Sized>(&mut self, problem: &P, stats: &mut Stats, t: f64) {
        let delta = f64::EPSILON.sqrt() * t.abs().max(1.0);
        problem.rhs(t, &self.y, &mut self.f0);
        problem.rhs(t + delta, &self.y, &mut self.ft);
        stats.rhs_evals += 2;
        for i in 0..self.ft.len() {
            self.ft[i] = (self.ft[i] - self.f0[i]) / delta;
        }
    }
}

impl<P: Problem + ?Sized> Stepper<P> for RosenbrockStepper {
    fn control_order(&self) -> usize {
        self.control_order
    }

    fn state(&self) -> &[f64] {
        &self.y
    }

    fn proposed(&self) -> &[f64] {
        &self.y_new
    }

    fn start(&mut self, _problem: &P, _stats: &mut Stats, _t: f64, y: &[f64]) {
        self.y.copy_from_slice(y);
        self.y_new.copy_from_slice(y);
        self.linear.invalidate();
        self.last_h = 0.0;
    }

    fn attempt(&mut self, problem: &P, stats: &mut Stats, t: f64, h: f64) -> StepOutcome {
        let stages = self.tableau.stages;
        let dim = self.y.len();

        if self.needs_time_derivative {
            self.time_derivative(problem, stats, t);
        }

        for i in 0..stages {
            let scale = h * self.diagonal[i];
            // The age limit is `never`, because the Jacobian is invalidated by
            // `commit` instead. Within one step, retries included, the point it
            // belongs to has not moved, so one evaluation serves them all.
            if !self
                .linear
                .prepare(problem, stats, t, &self.y, scale, u32::MAX)
            {
                return StepOutcome::failed();
            }

            self.argument.copy_from_slice(&self.y);
            for j in 0..i {
                let coefficient = self.tableau.a[(i, j)];
                if coefficient != 0.0 {
                    simd::axpy(coefficient, &self.s[j], &mut self.argument);
                }
            }

            problem.rhs(t + self.tableau.c[i] * h, &self.argument, &mut self.rhs);
            stats.rhs_evals += 1;

            for j in 0..i {
                let coefficient = self.tableau.c_matrix[(i, j)];
                if coefficient != 0.0 {
                    simd::axpy(coefficient / h, &self.s[j], &mut self.rhs);
                }
            }
            if self.needs_time_derivative {
                simd::axpy(h * self.tableau.d[i], &self.ft, &mut self.rhs);
            }

            // The system as written has `I/(h gamma)` on the diagonal. Scaling
            // it up gives `I - h gamma J`, the matrix the factorization is for.
            simd::scale(scale, &mut self.rhs);
            if !self.linear.solve(&mut self.rhs) {
                return StepOutcome::failed();
            }
            stats.linear_solves += 1;
            self.s[i].copy_from_slice(&self.rhs);
        }

        simd::combine(&self.y, 1.0, &self.tableau.m[..stages], &self.s, &mut self.y_new);
        if !self.y_new.iter().all(|v| v.is_finite()) {
            return StepOutcome::failed();
        }

        let error = match &self.tableau.e {
            Some(e) => {
                simd::combine_into(1.0, &e[..stages], &self.s, &mut self.err);
                simd::error_scale(self.atol, self.rtol, &self.y, &self.y_new, &mut self.scale);
                simd::weighted_rms(&self.err[..dim], &self.scale[..dim])
            }
            None => 0.0,
        };

        self.last_h = h;
        StepOutcome { ok: true, error }
    }

    /// The step is taken, so the Jacobian belongs to a point that no longer
    /// exists.
    ///
    /// An implicit Runge-Kutta method may keep a stale Jacobian for a while:
    /// there it only steers a Newton iteration, and a worse steer costs
    /// iterations, not accuracy. Here it is the method. Replacing `f'(y_n)` by
    /// something else does not solve the same equations less well, it defines a
    /// different scheme, a W-method, whose order is generally lower than the one
    /// the coefficients were built for. So it is thrown away every step, which
    /// is the price a ROW method charges and the reason it never iterates.
    fn commit(&mut self, _t: f64, _h: f64) {
        self.y.copy_from_slice(&self.y_new);
        self.linear.invalidate();
    }

    fn reject(&mut self, _h: f64) {}

    /// A Rosenbrock step leaves no derivative at its own end point, so there is
    /// nothing to interpolate with that would beat the straight line the driver
    /// falls back to. Published dense output exists for some of these methods
    /// and would go here.
    fn interpolate(&self, _theta: f64, _out: &mut [f64]) -> bool {
        false
    }

    fn is_adaptive(&self) -> bool {
        self.tableau.e.is_some()
    }
}
