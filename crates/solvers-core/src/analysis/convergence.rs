//! Empirical convergence order.
//!
//! The order conditions say what the order *should* be. This measures what it
//! actually is: fixed step runs over a sequence of step sizes, compared against
//! a closed form solution where one exists and against a tightly integrated
//! reference otherwise. The two numbers disagreeing is informative, most often
//! because of order reduction on a stiff problem.

use crate::method::Method;
use crate::ode::{self, Options};
use crate::problem::{Problem, Stats};
use crate::problems::TestProblem;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct ConvergencePoint {
    pub h: f64,
    /// Relative error in the maximum norm at the end of the interval.
    pub error: f64,
    pub rhs_evals: u64,
    pub steps: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConvergenceStudy {
    pub method: String,
    pub problem: String,
    pub points: Vec<ConvergencePoint>,
    /// Slope of the fit through the usable points.
    pub estimated_order: f64,
    /// Slope between the last two usable points, which is closer to the
    /// asymptotic rate when the coarse steps are still preasymptotic.
    pub local_order: f64,
    /// True when the comparison used a closed form solution.
    pub exact_reference: bool,
}

/// Relative error in the maximum norm.
fn relative_error(computed: &[f64], reference: &[f64]) -> f64 {
    let mut worst = 0.0f64;
    for i in 0..reference.len() {
        let scale = reference[i].abs().max(1e-12);
        worst = worst.max((computed[i] - reference[i]).abs() / scale);
    }
    worst
}

/// Integrate a problem to high accuracy to obtain a reference end state.
pub fn reference_solution<P: Problem + ?Sized>(
    method: &Method,
    problem: &P,
    t_span: (f64, f64),
    y0: &[f64],
) -> Option<Vec<f64>> {
    let mut options = Options::with_tolerances(1e-13, 1e-15);
    options.max_steps = 10_000_000;
    let solution = ode::integrate(method, problem, t_span, y0, &options);
    if solution.succeeded() {
        solution.last().cloned()
    } else {
        None
    }
}

/// Run a fixed step convergence study.
///
/// `reference` is used only when the problem has no closed form solution.
pub fn study(
    method: &Method,
    problem: &dyn TestProblem,
    step_sizes: &[f64],
    reference_method: Option<&Method>,
) -> ConvergenceStudy {
    let t_span = problem.t_span();
    let y0 = problem.y0();

    let (target, exact_reference) = match problem.exact(t_span.1) {
        Some(exact) => (Some(exact), true),
        None => {
            let reference = reference_method
                .and_then(|m| reference_solution(m, problem, t_span, &y0));
            (reference, false)
        }
    };

    let mut points = Vec::new();
    if let Some(target) = &target {
        for &h in step_sizes {
            let mut options = Options::fixed_step(h);
            options.max_steps = 20_000_000;
            // The implicit solves must be far more accurate than the
            // discretization error, otherwise the measurement reports the
            // iteration tolerance instead of the order of the method.
            options.rtol = 1e-12;
            options.atol = 1e-14;
            options.nonlinear.tolerance = 1e-2;
            options.nonlinear.max_iterations = 50;
            options.max_jacobian_age = 0;
            let solution = ode::integrate(method, problem, t_span, &y0, &options);
            let error = match solution.last() {
                Some(y) if solution.succeeded() => relative_error(y, target),
                _ => f64::INFINITY,
            };
            points.push(ConvergencePoint {
                h,
                error,
                rhs_evals: solution.stats.rhs_evals,
                steps: solution.stats.accepted,
            });
        }
    }

    let (estimated_order, local_order) = fit_order(&points);
    ConvergenceStudy {
        method: method.id.clone(),
        problem: problem.id().to_string(),
        points,
        estimated_order,
        local_order,
        exact_reference,
    }
}

/// Least squares slope of `log(error)` against `log(h)`.
///
/// Points that hit round off or diverged are dropped, because including them
/// bends the fit and hides the real rate.
pub fn fit_order(points: &[ConvergencePoint]) -> (f64, f64) {
    let usable: Vec<&ConvergencePoint> = points
        .iter()
        .filter(|p| p.error.is_finite() && p.error > 5e-13 && p.error < 0.5)
        .collect();
    if usable.len() < 2 {
        return (f64::NAN, f64::NAN);
    }

    let n = usable.len() as f64;
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    for point in &usable {
        let x = point.h.ln();
        let y = point.error.ln();
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denominator = n * sxx - sx * sx;
    let slope = if denominator.abs() < 1e-30 {
        f64::NAN
    } else {
        (n * sxy - sx * sy) / denominator
    };

    let last = usable[usable.len() - 1];
    let second_last = usable[usable.len() - 2];
    let local = (last.error.ln() - second_last.error.ln()) / (last.h.ln() - second_last.h.ln());

    (slope, local)
}

/// Default step size ladder: halving from `coarse` for `count` levels.
pub fn step_ladder(coarse: f64, count: usize) -> Vec<f64> {
    (0..count).map(|i| coarse / 2f64.powi(i as i32)).collect()
}

/// Total work of a run, for the record.
pub fn work(stats: &Stats) -> u64 {
    stats.rhs_evals
}
