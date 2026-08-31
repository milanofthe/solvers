//! The additive (IMEX) stepper.
//!
//! A step is one pass over the shared stages. At stage `i` everything already
//! computed is collected into a base point,
//!
//! ```text
//! base_i = y_n + h sum_{j<i} ( aE_ij kE_j + aI_ij kI_j )
//! ```
//!
//! and the only unknown left is the implicit half's own stage derivative, since
//! the explicit half contributes nothing on its diagonal:
//!
//! ```text
//! kI_i = f_I(t_n + cI_i h, base_i + h aI_ii kI_i)
//! Y_i  = base_i + h aI_ii kI_i
//! kE_i = f_E(t_n + cE_i h, Y_i)
//! ```
//!
//! The two abscissae are read separately because the halves are allowed to
//! differ on them, and some published pairs do.
//!
//! That first equation is the stage equation of a diagonally implicit method
//! with `f_I` in place of `f`, so it is solved by the same residual the DIRK
//! path uses, handed a view of the problem that answers with its implicit half.
//! One factorization serves every stage wherever the implicit half has a
//! constant diagonal, which every published pair does.
//!
//! A problem that states no splitting gives the whole right hand side to the
//! implicit half. The pair then runs as its implicit tableau alone, which is
//! what an IMEX method degenerates to when there is nothing to be explicit
//! about, and is the right answer rather than a special case.

use super::newton_matrix::NewtonMatrix;
use super::rk::DirkResidual;
use super::{Options, StepOutcome, Stepper};
use crate::method::{AdditiveTableau, RkRuntime};
use crate::nonlinear::NonlinearSolver;
use crate::problem::{Problem, Stats};
use crate::simd;

/// The implicit half of a problem, seen as a problem in its own right.
///
/// The Newton solve at a stage is over `f_I` alone, so it needs a right hand
/// side and a Jacobian that are `f_I`'s. Wrapping rather than branching keeps
/// the stage solve identical to the one the diagonally implicit methods use.
struct ImplicitPart<'a, P: Problem + ?Sized>(&'a P);

impl<'a, P: Problem + ?Sized> Problem for ImplicitPart<'a, P> {
    fn dim(&self) -> usize {
        self.0.dim()
    }
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        self.0.rhs_implicit(t, y, dy);
    }
}

pub struct AdditiveStepper {
    explicit: RkRuntime,
    implicit: RkRuntime,
    control_order: usize,
    dim: usize,
    rtol: f64,
    atol: f64,
    max_jacobian_age: u32,

    y: Vec<f64>,
    y_new: Vec<f64>,
    /// Stage derivatives of the two halves.
    ke: Vec<Vec<f64>>,
    ki: Vec<Vec<f64>>,
    base: Vec<f64>,
    stage_value: Vec<f64>,
    err: Vec<f64>,
    scale: Vec<f64>,
    work: Vec<f64>,
    stage: Vec<f64>,
    guess: Vec<f64>,

    linear: NewtonMatrix,
    nonlinear: NonlinearSolver,
    last_h: f64,
}

impl AdditiveStepper {
    pub fn new(
        pair: &AdditiveTableau,
        declared_order: usize,
        embedded_order: Option<usize>,
        dim: usize,
        options: &Options,
    ) -> AdditiveStepper {
        let explicit = pair.explicit.runtime();
        let implicit = pair.implicit.runtime();
        let s = explicit.stages;
        let control_order = match embedded_order {
            Some(m) => declared_order.min(m).max(1),
            None => declared_order.max(1),
        };
        AdditiveStepper {
            explicit,
            implicit,
            control_order,
            dim,
            rtol: options.rtol,
            atol: options.atol,
            max_jacobian_age: options.max_jacobian_age,
            y: vec![0.0; dim],
            y_new: vec![0.0; dim],
            ke: vec![vec![0.0; dim]; s],
            ki: vec![vec![0.0; dim]; s],
            base: vec![0.0; dim],
            stage_value: vec![0.0; dim],
            err: vec![0.0; dim],
            scale: vec![1.0; dim],
            work: vec![0.0; dim],
            stage: vec![0.0; dim],
            guess: vec![0.0; dim],
            linear: NewtonMatrix::new(dim),
            nonlinear: NonlinearSolver::new(options.nonlinear),
            last_h: 0.0,
        }
    }

    /// `y_n + h sum_{j<i} ( aE_ij kE_j + aI_ij kI_j )`.
    fn base_point(&mut self, i: usize, h: f64) {
        self.base.copy_from_slice(&self.y);
        for j in 0..i {
            let ae = self.explicit.a[(i, j)] * h;
            let ai = self.implicit.a[(i, j)] * h;
            if ae != 0.0 {
                for d in 0..self.dim {
                    self.base[d] += ae * self.ke[j][d];
                }
            }
            if ai != 0.0 {
                for d in 0..self.dim {
                    self.base[d] += ai * self.ki[j][d];
                }
            }
        }
    }

    /// The solution and the error estimate, from both halves.
    fn finalize(&mut self, h: f64) -> StepOutcome {
        let s = self.explicit.stages;
        simd::combine(&self.y, h, &self.explicit.b[..s], &self.ke, &mut self.y_new);
        for i in 0..s {
            let bi = self.implicit.b[i] * h;
            if bi != 0.0 {
                for d in 0..self.dim {
                    self.y_new[d] += bi * self.ki[i][d];
                }
            }
        }
        if !self.y_new.iter().all(|v| v.is_finite()) {
            return StepOutcome::failed();
        }

        // Both halves carry their own difference to the embedded weights, and
        // the estimate is the sum. It goes through the same filter the stiff
        // one step methods use, or it would report the stiff modes the implicit
        // half is damping perfectly as error.
        let error = match (&self.explicit.e, &self.implicit.e) {
            (Some(ee), Some(ei)) => {
                simd::combine_into(h, &ee[..s], &self.ke, &mut self.err);
                for i in 0..s {
                    let e = ei[i] * h;
                    if e != 0.0 {
                        for d in 0..self.dim {
                            self.err[d] += e * self.ki[i][d];
                        }
                    }
                }
                self.linear.solve(&mut self.err);
                simd::error_scale(self.atol, self.rtol, &self.y, &self.y_new, &mut self.scale);
                simd::weighted_rms(&self.err, &self.scale)
            }
            _ => 0.0,
        };
        self.last_h = h;
        StepOutcome { ok: true, error }
    }
}

impl<P: Problem + ?Sized> Stepper<P> for AdditiveStepper {
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
        self.last_h = 0.0;
    }

    fn attempt(&mut self, problem: &P, stats: &mut Stats, t: f64, h: f64) -> StepOutcome {
        let s = self.explicit.stages;
        let implicit_view = ImplicitPart(problem);

        if let Some(gamma) = self.implicit.gamma {
            if !self.linear.prepare(
                &implicit_view,
                stats,
                t,
                &self.y,
                h * gamma,
                self.max_jacobian_age,
            ) {
                return StepOutcome::failed();
            }
        }

        stats.rhs_evals += 1;
        problem.rhs_implicit(t, &self.y, &mut self.guess);
        simd::error_scale(self.atol, self.rtol, &self.y, &self.y, &mut self.scale);
        let full_newton = self.nonlinear.config.kind == crate::nonlinear::SolverKind::Newton;
        let max_age = self.max_jacobian_age;

        for i in 0..s {
            self.base_point(i, h);
            let ti = t + self.implicit.c[i] * h;
            let diagonal = self.implicit.a[(i, i)];

            if diagonal.abs() < 1e-14 {
                self.stage_value.copy_from_slice(&self.base);
                stats.rhs_evals += 1;
                let mut out = std::mem::take(&mut self.ki[i]);
                problem.rhs_implicit(ti, &self.stage_value, &mut out);
                self.ki[i] = out;
            } else {
                let mut residual = DirkResidual {
                    problem: &implicit_view,
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
                self.ki[i].copy_from_slice(&stage);
                self.guess.copy_from_slice(&stage);
                self.stage = stage;

                for d in 0..self.dim {
                    self.stage_value[d] = self.base[d] + h * diagonal * self.ki[i][d];
                }
            }

            stats.rhs_evals += 1;
            let mut out = std::mem::take(&mut self.ke[i]);
            problem.rhs_explicit(t + self.explicit.c[i] * h, &self.stage_value, &mut out);
            self.ke[i] = out;

            if !self.ki[i].iter().all(|v| v.is_finite())
                || !self.ke[i].iter().all(|v| v.is_finite())
            {
                return StepOutcome::failed();
            }
        }
        self.finalize(h)
    }

    fn commit(&mut self, _t: f64, _h: f64) {
        self.y.copy_from_slice(&self.y_new);
    }

    fn reject(&mut self, _h: f64) {}

    fn interpolate(&self, _theta: f64, _out: &mut [f64]) -> bool {
        false
    }

    fn is_adaptive(&self) -> bool {
        self.explicit.e.is_some() && self.implicit.e.is_some()
    }
}
