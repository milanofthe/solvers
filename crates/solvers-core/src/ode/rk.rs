//! The Runge-Kutta stepper.
//!
//! One stepper serves every tableau. The sparsity pattern of `A` decides how a
//! step is computed:
//!
//! * strictly lower triangular: the stages evaluate in sequence,
//! * lower triangular: one implicit solve of size `n` per stage, sharing a
//!   single factorization when the diagonal is constant (SDIRK, ESDIRK),
//! * dense: one coupled implicit solve of size `s * n` (Gauss, Radau, Lobatto).
//!
//! The coupled case is where the vectorized kernels earn their keep: the stage
//! coupling is a small dense matrix acting on `s` long vectors, which
//! `simd::stage_transform` evaluates blocked over components.

use super::decouple::{DecoupledLinear, StageDecoupling};
use super::newton_matrix::NewtonMatrix;
use super::{Options, StepOutcome, Stepper};
use crate::linalg::{Lu, Matrix};
use crate::method::{RkRuntime, RkTableau, Structure};
use crate::nonlinear::{NonlinearSolver, Residual};
use crate::problem::{Problem, Stats};
use crate::simd;

pub struct RkStepper {
    tableau: RkRuntime,
    control_order: usize,
    dim: usize,
    rtol: f64,
    atol: f64,
    max_jacobian_age: u32,

    y: Vec<f64>,
    y_new: Vec<f64>,
    /// Stage derivatives.
    k: Vec<Vec<f64>>,
    base: Vec<f64>,
    err: Vec<f64>,
    scale: Vec<f64>,
    /// Scratch buffers for stage values and implicit stage unknowns.
    work: Vec<f64>,
    stage: Vec<f64>,
    guess: Vec<f64>,
    stage_scale: Vec<f64>,
    stage_f: Vec<Vec<f64>>,
    stage_tmp: Vec<Vec<f64>>,
    /// `f(t, y)` at the start of the step when it is already known.
    f_start: Option<Vec<f64>>,
    /// True when stage one is `f(t, y)`, which makes FSAL reuse legal.
    trivial_first_stage: bool,
    /// Whether the embedded error estimate is filtered through the iteration
    /// matrix. Only meaningful once there is one, so only for implicit methods.
    filter_error: bool,
    /// Whether the error estimate comes from step doubling rather than from an
    /// embedded pair.
    richardson: bool,
    y_saved: Vec<f64>,
    y_big: Vec<f64>,

    linear: NewtonMatrix,
    nonlinear: NonlinearSolver,

    /// Coupled stage system for fully implicit tableaux.
    coupled: Option<CoupledLinear>,
    z: Vec<f64>,

    last_h: f64,
    have_step: bool,
}

/// Linear model of the coupled `s * n` stage system.
///
/// When the stage matrix diagonalizes, the system is solved as `s` decoupled
/// problems of size `n` instead; the dense factorization is the fallback for a
/// defective or singular `A`.
struct CoupledLinear {
    jacobian: Matrix<f64>,
    lu: Option<Lu<f64>>,
    decoupled: Option<DecoupledLinear>,
    jacobian_valid: bool,
    factored_h: f64,
    age: u32,
}

impl RkStepper {
    pub fn new(
        tableau: &RkTableau,
        declared_order: usize,
        embedded_order: Option<usize>,
        dim: usize,
        options: &Options,
    ) -> RkStepper {
        let runtime = tableau.runtime();
        let s = runtime.stages;

        let control_order = match embedded_order {
            Some(m) => declared_order.min(m).max(1),
            None => declared_order.max(1),
        };

        let trivial_first_stage = runtime.c[0].abs() < 1e-14
            && (0..s).all(|j| runtime.a[(0, j)].abs() < 1e-14);
        let filter_error = runtime.structure == Structure::DiagonallyImplicit;
        // Step doubling only makes sense when the driver is actually going to use
        // the estimate; a fixed step run would just pay three times the cost.
        let richardson = options.adaptive && options.richardson && runtime.e.is_none();

        let coupled = if runtime.structure == Structure::FullyImplicit {
            let decoupled = if options.decouple_stages {
                StageDecoupling::new(&runtime.a).map(|d| DecoupledLinear::new(d, dim))
            } else {
                None
            };
            Some(CoupledLinear {
                jacobian: Matrix::zeros(dim, dim),
                lu: None,
                decoupled,
                jacobian_valid: false,
                factored_h: f64::NAN,
                age: 0,
            })
        } else {
            None
        };

        RkStepper {
            tableau: runtime,
            control_order,
            dim,
            rtol: options.rtol,
            atol: options.atol,
            max_jacobian_age: options.max_jacobian_age,
            y: vec![0.0; dim],
            y_new: vec![0.0; dim],
            k: vec![vec![0.0; dim]; s],
            base: vec![0.0; dim],
            err: vec![0.0; dim],
            scale: vec![1.0; dim],
            work: vec![0.0; dim],
            stage: vec![0.0; dim],
            guess: vec![0.0; dim],
            stage_scale: vec![1.0; s * dim],
            stage_f: vec![vec![0.0; dim]; s],
            stage_tmp: vec![vec![0.0; dim]; s],
            f_start: None,
            trivial_first_stage,
            filter_error,
            richardson,
            y_saved: vec![0.0; dim],
            y_big: vec![0.0; dim],
            linear: NewtonMatrix::new(dim),
            nonlinear: NonlinearSolver::new(options.nonlinear),
            coupled,
            z: vec![0.0; s * dim],
            last_h: 0.0,
            have_step: false,
        }
    }

    pub fn tableau(&self) -> &RkRuntime {
        &self.tableau
    }

    /// Finish a step: form the solution and the error estimate.
    ///
    /// For an implicit method the raw embedded difference is passed through
    /// `(I - h*gamma*J)^{-1}` first. Without that filter the estimate inherits
    /// the stiff eigenvalues of the problem and reports an error of order
    /// `h*lambda` for modes the method is in fact damping perfectly, which
    /// throttles the step size to the explicit stability limit and throws away
    /// the entire point of using an L-stable method.
    ///
    /// Reference: Hairer and Wanner, "Solving ODEs II", IV.8, the error
    /// estimate of RADAU5.
    fn finalize(&mut self, h: f64) -> StepOutcome {
        let s = self.tableau.stages;
        simd::combine(&self.y, h, &self.tableau.b[..s], &self.k, &mut self.y_new);

        if !self.y_new.iter().all(|v| v.is_finite()) {
            return StepOutcome::failed();
        }

        let error = match &self.tableau.e {
            Some(e) => {
                simd::combine_into(h, &e[..s], &self.k, &mut self.err);
                if self.filter_error {
                    self.linear.solve(&mut self.err);
                }
                simd::error_scale(self.atol, self.rtol, &self.y, &self.y_new, &mut self.scale);
                simd::weighted_rms(&self.err, &self.scale)
            }
            None => 0.0,
        };
        self.last_h = h;
        StepOutcome { ok: true, error }
    }

    /// One step with the method as written.
    fn attempt_once<P: Problem + ?Sized>(
        &mut self,
        problem: &P,
        stats: &mut Stats,
        t: f64,
        h: f64,
    ) -> StepOutcome {
        match self.tableau.structure {
            Structure::Explicit => self.attempt_explicit(problem, stats, t, h),
            Structure::DiagonallyImplicit => self.attempt_dirk(problem, stats, t, h),
            Structure::FullyImplicit => self.attempt_coupled(problem, stats, t, h),
        }
    }

    /// One step of `h` against two of `h/2`, which gives an error estimate for
    /// a method that has no embedded pair of its own.
    ///
    /// The difference between the two results is `(2^p - 1)` times the leading
    /// error term of the single step, so dividing by that yields an estimate of
    /// the same quality an embedded pair would give. It costs three steps
    /// instead of one, which is the price of making Gauss, Radau and Lobatto
    /// adaptive at all; a method that ships an embedded pair never takes this
    /// path.
    fn attempt_with_step_doubling<P: Problem + ?Sized>(
        &mut self,
        problem: &P,
        stats: &mut Stats,
        t: f64,
        h: f64,
    ) -> StepOutcome {
        let half = 0.5 * h;
        self.y_saved.copy_from_slice(&self.y);

        let outcome = self.attempt_once(problem, stats, t, h);
        if !outcome.ok {
            return outcome;
        }
        self.y_big.copy_from_slice(&self.y_new);

        // Two half steps from the same starting point. The cached start
        // derivative belongs to a different step size, so it is dropped.
        self.f_start = None;
        let first = self.attempt_once(problem, stats, t, half);
        if !first.ok {
            self.y.copy_from_slice(&self.y_saved);
            self.f_start = None;
            return StepOutcome::failed();
        }
        self.y.copy_from_slice(&self.y_new);
        self.f_start = None;
        let second = self.attempt_once(problem, stats, t + half, half);
        self.y.copy_from_slice(&self.y_saved);
        self.f_start = None;
        if !second.ok {
            return StepOutcome::failed();
        }

        let order = self.control_order.max(1) as i32;
        let denominator = 2f64.powi(order) - 1.0;
        for i in 0..self.dim {
            self.err[i] = (self.y_new[i] - self.y_big[i]) / denominator;
        }
        simd::error_scale(self.atol, self.rtol, &self.y, &self.y_new, &mut self.scale);
        let error = simd::weighted_rms(&self.err, &self.scale);
        self.last_h = h;
        StepOutcome { ok: true, error }
    }

    fn attempt_explicit<P: Problem + ?Sized>(
        &mut self,
        problem: &P,
        stats: &mut Stats,
        t: f64,
        h: f64,
    ) -> StepOutcome {
        let s = self.tableau.stages;
        for i in 0..s {
            if i == 0 && self.trivial_first_stage {
                if let Some(f0) = &self.f_start {
                    self.k[0].copy_from_slice(f0);
                    continue;
                }
            }
            let row: Vec<f64> = (0..i).map(|j| self.tableau.a[(i, j)]).collect();
            simd::combine(&self.y, h, &row, &self.k, &mut self.base);
            stats.rhs_evals += 1;
            let ti = t + self.tableau.c[i] * h;
            let mut out = std::mem::take(&mut self.k[i]);
            problem.rhs(ti, &self.base, &mut out);
            self.k[i] = out;
            if !self.k[i].iter().all(|v| v.is_finite()) {
                return StepOutcome::failed();
            }
        }
        self.finalize(h)
    }

    fn attempt_dirk<P: Problem + ?Sized>(
        &mut self,
        problem: &P,
        stats: &mut Stats,
        t: f64,
        h: f64,
    ) -> StepOutcome {
        let s = self.tableau.stages;

        // A constant diagonal means one factorization serves every stage, which
        // is the reason SDIRK and ESDIRK methods exist at all.
        if let Some(gamma) = self.tableau.gamma {
            if !self
                .linear
                .prepare(problem, stats, t, &self.y, h * gamma, self.max_jacobian_age)
            {
                return StepOutcome::failed();
            }
        }

        if self.f_start.is_none() {
            let mut f0 = vec![0.0; self.dim];
            stats.rhs_evals += 1;
            problem.rhs(t, &self.y, &mut f0);
            self.f_start = Some(f0);
        }
        self.guess.copy_from_slice(self.f_start.as_ref().unwrap());

        simd::error_scale(self.atol, self.rtol, &self.y, &self.y, &mut self.scale);
        let full_newton = self.nonlinear.config.kind == crate::nonlinear::SolverKind::Newton;
        let max_age = self.max_jacobian_age;

        for i in 0..s {
            let row: Vec<f64> = (0..i).map(|j| self.tableau.a[(i, j)]).collect();
            simd::combine(&self.y, h, &row, &self.k, &mut self.base);
            let ti = t + self.tableau.c[i] * h;
            let diagonal = self.tableau.a[(i, i)];

            if diagonal.abs() < 1e-14 {
                stats.rhs_evals += 1;
                let mut out = std::mem::take(&mut self.k[i]);
                problem.rhs(ti, &self.base, &mut out);
                self.k[i] = out;
            } else {
                // The buffers move into the residual and back, so the nonlinear
                // solver and the linear model can be borrowed at the same time.
                let mut residual = DirkResidual {
                    problem,
                    stats,
                    linear: &mut self.linear,
                    anchor_t: t,
                    anchor_y: std::mem::take(&mut self.y),
                    t: ti,
                    base: std::mem::take(&mut self.base),
                    h_gamma: h * diagonal,
                    max_age,
                    work: std::mem::take(&mut self.work),
                    full_newton,
                };

                let mut stage = std::mem::take(&mut self.stage);
                stage.copy_from_slice(&self.guess);
                self.nonlinear.reset();
                let outcome = self.nonlinear.solve(&mut residual, &mut stage, &self.scale);

                self.y = residual.anchor_y;
                self.base = residual.base;
                self.work = residual.work;
                outcome.record(stats);

                if !outcome.converged() {
                    self.stage = stage;
                    return StepOutcome::failed();
                }
                self.k[i].copy_from_slice(&stage);
                self.guess.copy_from_slice(&stage);
                self.stage = stage;
            }

            if !self.k[i].iter().all(|v| v.is_finite()) {
                return StepOutcome::failed();
            }
        }
        self.finalize(h)
    }

    fn attempt_coupled<P: Problem + ?Sized>(
        &mut self,
        problem: &P,
        stats: &mut Stats,
        t: f64,
        h: f64,
    ) -> StepOutcome {
        let s = self.tableau.stages;
        let n = self.dim;

        // Predict the stage increments from the previous stage derivatives. For
        // a smooth solution this beats starting at zero and saves an iteration
        // or two per step.
        if self.have_step {
            for i in 0..s {
                for p in 0..n {
                    self.z[i * n + p] = h * self.tableau.c[i] * self.k[i][p];
                }
            }
        } else {
            for v in self.z.iter_mut() {
                *v = 0.0;
            }
        }

        simd::error_scale(self.atol, self.rtol, &self.y, &self.y, &mut self.scale);
        for idx in 0..s * n {
            self.stage_scale[idx] = self.scale[idx % n];
        }

        let mut residual = CoupledResidual {
            problem,
            stats,
            linear: self.coupled.as_mut().expect("coupled linear model missing"),
            tableau: &self.tableau,
            y: std::mem::take(&mut self.y),
            t,
            h,
            max_age: self.max_jacobian_age,
            f: std::mem::take(&mut self.stage_f),
            tmp: std::mem::take(&mut self.stage_tmp),
            work: std::mem::take(&mut self.work),
        };

        self.nonlinear.reset();
        let outcome = self.nonlinear.solve(&mut residual, &mut self.z, &self.stage_scale);

        self.y = residual.y;
        self.stage_f = residual.f;
        self.stage_tmp = residual.tmp;
        self.work = residual.work;
        outcome.record(stats);

        if !outcome.converged() {
            return StepOutcome::failed();
        }
        for i in 0..s {
            self.k[i].copy_from_slice(&self.stage_f[i]);
        }
        self.finalize(h)
    }
}

impl<P: Problem + ?Sized> Stepper<P> for RkStepper {
    fn control_order(&self) -> usize {
        self.control_order
    }

    fn state(&self) -> &[f64] {
        &self.y
    }

    fn proposed(&self) -> &[f64] {
        &self.y_new
    }

    fn is_adaptive(&self) -> bool {
        self.tableau.e.is_some() || self.richardson
    }

    fn start(&mut self, _problem: &P, _stats: &mut Stats, _t: f64, y: &[f64]) {
        self.y.copy_from_slice(y);
        self.y_new.copy_from_slice(y);
        self.f_start = None;
        self.have_step = false;
        self.linear.invalidate();
        self.nonlinear.reset();
        if let Some(c) = &mut self.coupled {
            c.jacobian_valid = false;
            c.lu = None;
            if let Some(decoupled) = &mut c.decoupled {
                decoupled.invalidate();
            }
        }
    }

    fn attempt(&mut self, problem: &P, stats: &mut Stats, t: f64, h: f64) -> StepOutcome {
        if self.richardson {
            return self.attempt_with_step_doubling(problem, stats, t, h);
        }
        self.attempt_once(problem, stats, t, h)
    }

    fn commit(&mut self, _t: f64, _h: f64) {
        self.y.copy_from_slice(&self.y_new);
        self.have_step = true;
        if self.richardson {
            // The cached derivative belongs to the last half step.
            self.f_start = None;
        }
        self.linear.advance_age();
        if let Some(c) = &mut self.coupled {
            c.age += 1;
        }
        // The last stage of a stiffly accurate method is f at the new point.
        if self.tableau.stiffly_accurate {
            let last = self.tableau.stages - 1;
            match &mut self.f_start {
                Some(buffer) => buffer.copy_from_slice(&self.k[last]),
                None => self.f_start = Some(self.k[last].clone()),
            }
        } else {
            self.f_start = None;
        }
    }

    fn reject(&mut self, _h: f64) {
        // The stage values are stale but the Jacobian is still a fair model.
        self.nonlinear.reset();
    }

    fn interpolate(&self, theta: f64, out: &mut [f64]) -> bool {
        if self.richardson {
            // The stage values belong to the last half step, not to the whole
            // step, so they cannot be used to interpolate across it.
            return false;
        }
        let Some(dense) = &self.tableau.dense else {
            // Without a published interpolant, a cubic Hermite through the two
            // endpoint derivatives is available for free whenever the first
            // stage is `f(t, y)` and the last stage is `f(t + h, y_new)`, which
            // covers every FSAL and every stiffly accurate method.
            if self.trivial_first_stage && self.tableau.stiffly_accurate {
                let f0 = &self.k[0];
                let f1 = &self.k[self.tableau.stages - 1];
                let t2 = theta * theta;
                let t3 = t2 * theta;
                let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
                let h10 = t3 - 2.0 * t2 + theta;
                let h01 = -2.0 * t3 + 3.0 * t2;
                let h11 = t3 - t2;
                for i in 0..out.len() {
                    out[i] = h00 * self.y[i]
                        + h10 * self.last_h * f0[i]
                        + h01 * self.y_new[i]
                        + h11 * self.last_h * f1[i];
                }
                return true;
            }
            return false;
        };
        let s = self.tableau.stages;
        let degree = dense.cols();
        let mut weights = vec![0.0; s];
        for i in 0..s {
            let mut acc = 0.0;
            let mut power = theta;
            for j in 0..degree {
                acc += dense[(i, j)] * power;
                power *= theta;
            }
            weights[i] = acc;
        }
        simd::combine(&self.y, self.last_h, &weights, &self.k, out);
        true
    }
}

/// Stage equation of a diagonally implicit method,
/// `G(k) = k - f(t_i, base + h*a_ii*k)`.
struct DirkResidual<'a, P: Problem + ?Sized> {
    problem: &'a P,
    stats: &'a mut Stats,
    linear: &'a mut NewtonMatrix,
    /// Point the Jacobian is anchored at, the start of the step.
    anchor_t: f64,
    anchor_y: Vec<f64>,
    t: f64,
    base: Vec<f64>,
    h_gamma: f64,
    max_age: u32,
    work: Vec<f64>,
    full_newton: bool,
}

impl<'a, P: Problem + ?Sized> Residual for DirkResidual<'a, P> {
    fn dim(&self) -> usize {
        self.base.len()
    }

    fn eval(&mut self, z: &[f64], r: &mut [f64]) {
        for i in 0..self.work.len() {
            self.work[i] = self.base[i] + self.h_gamma * z[i];
        }
        self.stats.rhs_evals += 1;
        self.problem.rhs(self.t, &self.work, r);
        for i in 0..r.len() {
            r[i] = z[i] - r[i];
        }
    }

    fn factor(&mut self, z: &[f64]) -> bool {
        if self.full_newton {
            for i in 0..self.work.len() {
                self.work[i] = self.base[i] + self.h_gamma * z[i];
            }
            self.linear.invalidate();
            let work = std::mem::take(&mut self.work);
            let ok = self
                .linear
                .prepare(self.problem, self.stats, self.t, &work, self.h_gamma, 0);
            self.work = work;
            return ok;
        }
        self.linear.prepare(
            self.problem,
            self.stats,
            self.anchor_t,
            &self.anchor_y,
            self.h_gamma,
            self.max_age,
        )
    }

    fn solve(&mut self, rhs: &mut [f64]) -> bool {
        self.linear.solve(rhs)
    }

    fn refresh(&mut self) {
        self.linear.invalidate();
    }
}

/// Coupled stage system of a fully implicit method.
///
/// The unknown is the block vector of stage increments `z_i = Y_i - y`, and
/// `G(Z)_i = z_i - h * sum_j a_ij f(t + c_j h, y + z_j)`.
struct CoupledResidual<'a, P: Problem + ?Sized> {
    problem: &'a P,
    stats: &'a mut Stats,
    linear: &'a mut CoupledLinear,
    tableau: &'a RkRuntime,
    y: Vec<f64>,
    t: f64,
    h: f64,
    max_age: u32,
    /// Stage derivatives, kept so the caller can reuse them for the solution.
    f: Vec<Vec<f64>>,
    tmp: Vec<Vec<f64>>,
    work: Vec<f64>,
}

impl<'a, P: Problem + ?Sized> Residual for CoupledResidual<'a, P> {
    fn dim(&self) -> usize {
        self.tableau.stages * self.y.len()
    }

    fn eval(&mut self, z: &[f64], r: &mut [f64]) {
        let s = self.tableau.stages;
        let n = self.y.len();
        for j in 0..s {
            for p in 0..n {
                self.work[p] = self.y[p] + z[j * n + p];
            }
            self.stats.rhs_evals += 1;
            self.problem
                .rhs(self.t + self.tableau.c[j] * self.h, &self.work, &mut self.f[j]);
        }
        // tmp_i = h * sum_j a_ij f_j, the only place the stage coupling enters.
        simd::stage_transform(
            self.tableau.a.as_slice(),
            s,
            self.h,
            &self.f,
            &mut self.tmp,
        );
        for i in 0..s {
            for p in 0..n {
                r[i * n + p] = z[i * n + p] - self.tmp[i][p];
            }
        }
    }

    fn factor(&mut self, _z: &[f64]) -> bool {
        let n = self.y.len();
        let s = self.tableau.stages;

        let refresh = !self.linear.jacobian_valid || self.linear.age >= self.max_age;
        if refresh {
            self.stats.jacobian_evals += 1;
            if !self.problem.has_analytic_jacobian() {
                self.stats.rhs_evals += n as u64 + 1;
            }
            self.problem
                .jacobian(self.t, &self.y, &mut self.linear.jacobian);
            self.linear.jacobian_valid = true;
            self.linear.age = 0;
            self.linear.lu = None;
        }
        // The decoupled path: one factorization of size n per real eigenvalue
        // of A^{-1} and one per conjugate pair.
        if let Some(decoupled) = &mut self.linear.decoupled {
            if decoupled.is_factored_for(self.h) {
                return true;
            }
            if decoupled.factor(&self.linear.jacobian, self.h) {
                self.stats.lu_decompositions += decoupled.factorization_count() as u64;
                self.linear.factored_h = self.h;
                return true;
            }
            return false;
        }

        if self.linear.lu.is_some() && self.linear.factored_h == self.h {
            return true;
        }

        // M = I - h * (A kron J)
        let mut m = Matrix::<f64>::zeros(s * n, s * n);
        for i in 0..s {
            for j in 0..s {
                let a = self.tableau.a[(i, j)];
                if a == 0.0 {
                    continue;
                }
                let factor = -self.h * a;
                for p in 0..n {
                    for q in 0..n {
                        m[(i * n + p, j * n + q)] = factor * self.linear.jacobian[(p, q)];
                    }
                }
            }
        }
        for d in 0..s * n {
            m[(d, d)] += 1.0;
        }

        let lu = Lu::factor(m);
        if lu.is_singular() {
            self.linear.lu = None;
            return false;
        }
        self.stats.lu_decompositions += 1;
        self.linear.factored_h = self.h;
        self.linear.lu = Some(lu);
        true
    }

    fn solve(&mut self, rhs: &mut [f64]) -> bool {
        let h = self.h;
        if let Some(decoupled) = &mut self.linear.decoupled {
            return decoupled.solve(rhs, h);
        }
        match &self.linear.lu {
            Some(lu) => lu.solve_in_place(rhs),
            None => false,
        }
    }

    fn refresh(&mut self) {
        self.linear.jacobian_valid = false;
        self.linear.lu = None;
        if let Some(decoupled) = &mut self.linear.decoupled {
            decoupled.invalidate();
        }
    }
}
