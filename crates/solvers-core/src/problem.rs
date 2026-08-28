//! The initial value problem interface.

use crate::linalg::Matrix;

/// An explicit first order system `y' = f(t, y)`.
pub trait Problem {
    fn dim(&self) -> usize;

    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]);

    /// Whether `jacobian` is analytic. Purely informational, the stepper works
    /// either way, but the cost model reports the two cases separately.
    fn has_analytic_jacobian(&self) -> bool {
        false
    }

    /// Jacobian `df/dy`. The default is a first order forward difference.
    fn jacobian(&self, t: f64, y: &[f64], j: &mut Matrix<f64>) {
        finite_difference_jacobian(self, t, y, j)
    }
}

/// Forward difference Jacobian with the usual square root of eps scaling.
pub fn finite_difference_jacobian<P: Problem + ?Sized>(
    problem: &P,
    t: f64,
    y: &[f64],
    j: &mut Matrix<f64>,
) {
    let n = problem.dim();
    let mut f0 = vec![0.0; n];
    let mut f1 = vec![0.0; n];
    let mut yp = y.to_vec();
    problem.rhs(t, y, &mut f0);

    let sqrt_eps = f64::EPSILON.sqrt();
    for col in 0..n {
        let delta = sqrt_eps * y[col].abs().max(1.0) * if y[col] < 0.0 { -1.0 } else { 1.0 };
        let saved = yp[col];
        yp[col] = saved + delta;
        // Use the actually representable increment.
        let step = yp[col] - saved;
        problem.rhs(t, &yp, &mut f1);
        yp[col] = saved;
        for row in 0..n {
            j[(row, col)] = (f1[row] - f0[row]) / step;
        }
    }
}

/// Counters for the work a run actually performed.
///
/// These are what the cost analysis plots against the achieved accuracy, so
/// every implementation increments them at exactly one place.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Stats {
    pub steps: u64,
    pub accepted: u64,
    pub rejected: u64,
    pub rhs_evals: u64,
    pub jacobian_evals: u64,
    pub lu_decompositions: u64,
    pub linear_solves: u64,
    pub nonlinear_iterations: u64,
    pub nonlinear_failures: u64,
}

impl Stats {
    pub fn reset(&mut self) {
        *self = Stats::default();
    }
}

/// Problem wrapper that keeps the evaluation counters honest.
pub struct Counted<'a, P: Problem + ?Sized> {
    pub problem: &'a P,
    pub stats: &'a mut Stats,
}

impl<'a, P: Problem + ?Sized> Counted<'a, P> {
    pub fn new(problem: &'a P, stats: &'a mut Stats) -> Self {
        Counted { problem, stats }
    }

    pub fn dim(&self) -> usize {
        self.problem.dim()
    }

    pub fn rhs(&mut self, t: f64, y: &[f64], dy: &mut [f64]) {
        self.stats.rhs_evals += 1;
        self.problem.rhs(t, y, dy);
    }

    pub fn jacobian(&mut self, t: f64, y: &[f64], j: &mut Matrix<f64>) {
        self.stats.jacobian_evals += 1;
        if !self.problem.has_analytic_jacobian() {
            // A finite difference Jacobian costs n extra right hand sides.
            self.stats.rhs_evals += self.problem.dim() as u64 + 1;
        }
        self.problem.jacobian(t, y, j);
    }
}

/// Adapter that turns a closure into a `Problem`.
pub struct OdeFn<F> {
    dim: usize,
    f: F,
}

impl<F> OdeFn<F>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    pub fn new(dim: usize, f: F) -> Self {
        OdeFn { dim, f }
    }
}

impl<F> Problem for OdeFn<F>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    fn dim(&self) -> usize {
        self.dim
    }
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        (self.f)(t, y, dy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_difference_matches_analytic() {
        // y' = [-2 y0 + y1, y0 * y1]
        let p = OdeFn::new(2, |_t, y, dy| {
            dy[0] = -2.0 * y[0] + y[1];
            dy[1] = y[0] * y[1];
        });
        let y = [1.5, -0.7];
        let mut j = Matrix::zeros(2, 2);
        p.jacobian(0.0, &y, &mut j);
        assert!((j[(0, 0)] + 2.0).abs() < 1e-6);
        assert!((j[(0, 1)] - 1.0).abs() < 1e-6);
        assert!((j[(1, 0)] - y[1]).abs() < 1e-6);
        assert!((j[(1, 1)] - y[0]).abs() < 1e-6);
    }
}
