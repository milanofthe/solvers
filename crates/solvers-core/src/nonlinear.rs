//! Nonlinear solvers for implicit update equations.
//!
//! Every implicit method in this crate poses its stage equations as a root
//! problem `G(z) = 0` and hands it to one of the solvers here, so the choice of
//! iteration is independent of the choice of integration method.
//!
//! The problem object owns the linear model, because the iteration matrix of a
//! Runge-Kutta stage (`I - h*gamma*J`) is known analytically and rebuilding it
//! by finite differences would throw that away. A problem that offers no linear
//! model can still be solved by the fixed point iterations.
//!
//! References
//! ----------
//! * H. F. Walker, P. Ni, "Anderson acceleration for fixed-point iterations",
//!   SIAM J. Numer. Anal. 49(4), 2011, doi:10.1137/10078356X
//! * E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II",
//!   2nd ed., Springer 1996, doi:10.1007/978-3-642-05221-7

use crate::linalg::{Lu, Matrix};
use crate::simd;
use serde::{Deserialize, Serialize};

/// A root problem `G(z) = 0` with an optional linear model.
pub trait Residual {
    fn dim(&self) -> usize;

    /// Evaluate `G(z)` into `r`.
    fn eval(&mut self, z: &[f64], r: &mut [f64]);

    /// Rebuild and factor the iteration matrix at `z`.
    ///
    /// Returning `false` means no linear model is available, which restricts
    /// the caller to fixed point style iterations.
    fn factor(&mut self, _z: &[f64]) -> bool {
        false
    }

    /// Overwrite `rhs` with the solution of `M dz = rhs` for the factored `M`.
    fn solve(&mut self, _rhs: &mut [f64]) -> bool {
        false
    }

    /// Discard any cached derivative information, so the next `factor` rebuilds
    /// it from scratch. Called when the iteration stops contracting.
    fn refresh(&mut self) {}
}

/// Which iteration to run.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SolverKind {
    /// `z <- z - G(z)`. No linear model needed, contracts only for small steps.
    FixedPoint,
    /// Anderson accelerated fixed point iteration.
    Anderson,
    /// Newton with the iteration matrix rebuilt every iteration.
    Newton,
    /// Newton with a frozen iteration matrix, refactored when convergence
    /// deteriorates. The standard choice for stiff integrators.
    ModifiedNewton,
    /// Modified Newton whose increments are Anderson accelerated.
    NewtonAnderson,
}

impl SolverKind {
    pub fn from_name(name: &str) -> Option<SolverKind> {
        Some(
            match name.trim().to_ascii_lowercase().replace(['-', '_', ' '], "").as_str() {
                "fixedpoint" | "fpi" | "picard" => SolverKind::FixedPoint,
                "anderson" => SolverKind::Anderson,
                "newton" | "fullnewton" => SolverKind::Newton,
                "modifiednewton" | "simplifiednewton" | "quasinewton" => SolverKind::ModifiedNewton,
                "newtonanderson" => SolverKind::NewtonAnderson,
                _ => return None,
            },
        )
    }

    pub fn name(self) -> &'static str {
        match self {
            SolverKind::FixedPoint => "fixed_point",
            SolverKind::Anderson => "anderson",
            SolverKind::Newton => "newton",
            SolverKind::ModifiedNewton => "modified_newton",
            SolverKind::NewtonAnderson => "newton_anderson",
        }
    }

    pub fn all() -> &'static [SolverKind] {
        use SolverKind::*;
        &[FixedPoint, Anderson, Newton, ModifiedNewton, NewtonAnderson]
    }

    pub fn needs_jacobian(self) -> bool {
        matches!(
            self,
            SolverKind::Newton | SolverKind::ModifiedNewton | SolverKind::NewtonAnderson
        )
    }
}

#[derive(Copy, Clone, Debug, Deserialize, Serialize)]
pub struct NonlinearConfig {
    pub kind: SolverKind,
    /// Convergence threshold on the scaled increment norm.
    pub tolerance: f64,
    pub max_iterations: u32,
    /// Anderson history depth.
    pub depth: usize,
    /// Refactor once the observed contraction rate exceeds this.
    pub max_contraction: f64,
    /// Tikhonov regularization of the Anderson least squares problem.
    pub regularization: f64,
}

impl Default for NonlinearConfig {
    fn default() -> Self {
        NonlinearConfig {
            kind: SolverKind::ModifiedNewton,
            tolerance: 1e-3,
            max_iterations: 20,
            depth: 5,
            max_contraction: 0.5,
            regularization: 1e-12,
        }
    }
}

/// Why an iteration stopped.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Status {
    Converged,
    /// Ran out of iterations while still contracting.
    MaxIterations,
    /// The contraction rate exceeded one.
    Diverged,
    /// The linear model could not be factored or applied.
    LinearFailure,
}

#[derive(Copy, Clone, Debug)]
pub struct Outcome {
    pub status: Status,
    pub iterations: u32,
    /// Back substitutions performed against the factored iteration matrix.
    pub linear_solves: u32,
    /// Scaled norm of the last increment.
    pub increment: f64,
    /// Observed contraction rate, `NaN` before two increments are available.
    pub rate: f64,
}

impl Outcome {
    /// Fold the work of this solve into the run statistics.
    pub fn record(&self, stats: &mut crate::problem::Stats) {
        stats.nonlinear_iterations += self.iterations as u64;
        stats.linear_solves += self.linear_solves as u64;
        if self.status != Status::Converged {
            stats.nonlinear_failures += 1;
        }
    }

    pub fn converged(&self) -> bool {
        self.status == Status::Converged
    }
}

/// Iteration state, reused across steps so the Anderson history and the frozen
/// factorization survive when they are still useful.
pub struct NonlinearSolver {
    pub config: NonlinearConfig,
    // Anderson history.
    x_hist: Vec<Vec<f64>>,
    f_hist: Vec<Vec<f64>>,
    scratch: Vec<f64>,
    increment: Vec<f64>,
}

impl NonlinearSolver {
    pub fn new(config: NonlinearConfig) -> NonlinearSolver {
        NonlinearSolver {
            config,
            x_hist: Vec::new(),
            f_hist: Vec::new(),
            scratch: Vec::new(),
            increment: Vec::new(),
        }
    }

    /// Forget the acceleration history.
    pub fn reset(&mut self) {
        self.x_hist.clear();
        self.f_hist.clear();
    }

    /// Drop only the acceleration history, keeping the factorization.
    pub fn reset_history(&mut self) {
        self.x_hist.clear();
        self.f_hist.clear();
    }

    /// Solve `G(z) = 0` starting from `z`.
    ///
    /// `scale` gives the per component error weights the increment norm is
    /// measured against, so the convergence test matches the accuracy request.
    pub fn solve<R: Residual + ?Sized>(
        &mut self,
        problem: &mut R,
        z: &mut [f64],
        scale: &[f64],
    ) -> Outcome {
        let n = problem.dim();
        self.scratch.resize(n, 0.0);
        self.increment.resize(n, 0.0);
        self.x_hist.clear();
        self.f_hist.clear();

        let kind = self.config.kind;
        let newtonish = kind.needs_jacobian();

        // `factor` is idempotent and cheap when nothing changed, so it is always
        // called; the residual decides whether real work is needed.
        if newtonish && !problem.factor(z) {
            return Outcome {
                status: Status::LinearFailure,
                iterations: 0,
                linear_solves: 0,
                increment: f64::INFINITY,
                rate: f64::NAN,
            };
        }

        let mut previous = f64::NAN;
        let mut rate = f64::NAN;
        let mut refreshed = false;
        let mut linear_solves = 0u32;

        for iteration in 1..=self.config.max_iterations {

            // Residual, then the raw increment dz with z_new = z - dz.
            let mut residual = vec![0.0; n];
            problem.eval(z, &mut residual);

            if newtonish {
                if kind == SolverKind::Newton && iteration > 1 {
                    if !problem.factor(z) {
                        return Outcome {
                            status: Status::LinearFailure,
                            iterations: iteration,
                            linear_solves,
                            increment: f64::INFINITY,
                            rate,
                        };
                    }
                }
                if !problem.solve(&mut residual) {
                    problem.refresh();
                    return Outcome {
                        status: Status::LinearFailure,
                        iterations: iteration,
                        linear_solves,
                        increment: f64::INFINITY,
                        rate,
                    };
                }
                linear_solves += 1;
            }
            // For the fixed point family the residual is already the increment.

            let accelerate = matches!(kind, SolverKind::Anderson | SolverKind::NewtonAnderson)
                && self.config.depth > 0;

            if accelerate {
                // Anderson works on the map g(z) = z - dz.
                for i in 0..n {
                    self.scratch[i] = z[i] - residual[i];
                }
                self.anderson(z, &self.scratch.clone(), &mut residual);
            }

            let magnitude = simd::weighted_rms(&residual, scale);
            for i in 0..n {
                z[i] -= residual[i];
            }

            if !magnitude.is_finite() {
                problem.refresh();
                return Outcome {
                    status: Status::Diverged,
                    iterations: iteration,
                    linear_solves,
                    increment: magnitude,
                    rate,
                };
            }

            if iteration > 1 && previous > 0.0 {
                rate = magnitude / previous;
                if rate >= 1.0 {
                    // Give a frozen Jacobian one chance to catch up before
                    // declaring the iteration lost.
                    if newtonish && !refreshed {
                        refreshed = true;
                        problem.refresh();
                        self.x_hist.clear();
                        self.f_hist.clear();
                        if problem.factor(z) {
                            previous = magnitude;
                            continue;
                        }
                    }
                    problem.refresh();
                    return Outcome {
                        status: Status::Diverged,
                        iterations: iteration,
                        linear_solves,
                        increment: magnitude,
                        rate,
                    };
                }
            }

            // Hairer's criterion: bound the remaining error of the iteration by
            // the geometric tail rather than by the last increment alone.
            let projected = if rate.is_finite() && rate < 1.0 {
                magnitude * rate / (1.0 - rate)
            } else {
                magnitude
            };
            if projected <= self.config.tolerance {
                if rate.is_finite() && rate > self.config.max_contraction {
                    // Converged, but slowly enough that the next step should
                    // start from a fresh factorization.
                    problem.refresh();
                }
                return Outcome {
                    status: Status::Converged,
                    iterations: iteration,
                    linear_solves,
                    increment: magnitude,
                    rate,
                };
            }

            previous = magnitude;
        }

        problem.refresh();
        Outcome {
            status: Status::MaxIterations,
            iterations: self.config.max_iterations,
            linear_solves,
            increment: previous,
            rate,
        }
    }

    /// Replace the plain increment by the Anderson accelerated one.
    ///
    /// `x` is the current iterate, `g` the plain fixed point image, and
    /// `increment` is overwritten so that the caller's `x - increment` becomes
    /// the accelerated iterate.
    fn anderson(&mut self, x: &[f64], g: &[f64], increment: &mut [f64]) {
        let n = x.len();
        let residual: Vec<f64> = (0..n).map(|i| g[i] - x[i]).collect();

        self.x_hist.push(g.to_vec());
        self.f_hist.push(residual.clone());
        while self.x_hist.len() > self.config.depth + 1 {
            self.x_hist.remove(0);
            self.f_hist.remove(0);
        }

        let m = self.f_hist.len().saturating_sub(1);
        if m == 0 {
            return;
        }

        // Least squares over the increments of the residual history.
        let mut normal = Matrix::<f64>::zeros(m, m);
        let mut rhs = vec![0.0; m];
        let mut df: Vec<Vec<f64>> = Vec::with_capacity(m);
        let mut dg: Vec<Vec<f64>> = Vec::with_capacity(m);
        for k in 0..m {
            df.push(
                (0..n)
                    .map(|i| self.f_hist[k + 1][i] - self.f_hist[k][i])
                    .collect(),
            );
            dg.push(
                (0..n)
                    .map(|i| self.x_hist[k + 1][i] - self.x_hist[k][i])
                    .collect(),
            );
        }
        for i in 0..m {
            for j in 0..m {
                normal[(i, j)] = simd::dot(&df[i], &df[j]);
            }
            normal[(i, i)] += self.config.regularization;
            rhs[i] = simd::dot(&df[i], &residual);
        }

        let Some(gamma) = Lu::factor(normal).solve(&rhs) else {
            return;
        };
        if gamma.iter().any(|v| !v.is_finite()) {
            return;
        }

        // x_accelerated = g - sum_k gamma_k * dg_k, so the increment relative to
        // the current iterate is x - x_accelerated.
        for i in 0..n {
            let mut correction = 0.0;
            for k in 0..m {
                correction += gamma[k] * dg[k][i];
            }
            increment[i] = x[i] - (g[i] - correction);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /// G(z) = z - cos(z), a contraction with a known root.
    struct Cosine {
        jac: f64,
    }

    impl Residual for Cosine {
        fn dim(&self) -> usize {
            1
        }
        fn eval(&mut self, z: &[f64], r: &mut [f64]) {
            r[0] = z[0] - z[0].cos();
        }
        fn factor(&mut self, z: &[f64]) -> bool {
            self.jac = 1.0 + z[0].sin();
            true
        }
        fn solve(&mut self, rhs: &mut [f64]) -> bool {
            if self.jac.abs() < 1e-300 {
                return false;
            }
            rhs[0] /= self.jac;
            true
        }
    }

    fn run(kind: SolverKind) -> (f64, u32) {
        let mut cfg = NonlinearConfig::default();
        cfg.kind = kind;
        cfg.tolerance = 1e-12;
        cfg.max_iterations = 200;
        let mut solver = NonlinearSolver::new(cfg);
        let mut problem = Cosine { jac: 1.0 };
        let mut z = [0.0];
        let outcome = solver.solve(&mut problem, &mut z, &[1.0]);
        assert!(outcome.converged(), "{kind:?} did not converge: {outcome:?}");
        (z[0], outcome.iterations)
    }

    #[test]
    fn every_solver_finds_the_dottie_number() {
        // The fixed point of cos, 0.739085133215...
        for kind in SolverKind::all() {
            let (root, _) = run(*kind);
            assert!((root - 0.739085133215160).abs() < 1e-9, "{kind:?} gave {root}");
        }
    }

    #[test]
    fn newton_beats_plain_fixed_point() {
        let (_, newton) = run(SolverKind::Newton);
        let (_, fixed) = run(SolverKind::FixedPoint);
        assert!(newton < fixed);
    }

    #[test]
    fn anderson_beats_plain_fixed_point() {
        let (_, anderson) = run(SolverKind::Anderson);
        let (_, fixed) = run(SolverKind::FixedPoint);
        assert!(anderson < fixed);
    }

    #[test]
    fn names_round_trip() {
        for kind in SolverKind::all() {
            assert_eq!(SolverKind::from_name(kind.name()), Some(*kind));
        }
    }
}
