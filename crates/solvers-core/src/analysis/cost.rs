//! Work precision.
//!
//! The question a work precision diagram answers is the only one that matters
//! when choosing a method: at the accuracy you actually need, which method gets
//! there for the least work. Each point is one adaptive run at one tolerance,
//! plotting achieved error against the work it took.
//!
//! Work is reported in its components rather than as a single number, because
//! the right weighting depends on the problem: a cheap right hand side with an
//! expensive Jacobian is a different regime from the reverse.

use crate::method::Method;
use crate::ode::{self, Options};
use crate::problem::{Problem, Stats};
use crate::problems::TestProblem;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct CostPoint {
    pub rtol: f64,
    pub atol: f64,
    /// Relative error in the maximum norm at the end of the interval.
    pub error: f64,
    pub rhs_evals: u64,
    pub jacobian_evals: u64,
    pub lu_decompositions: u64,
    pub linear_solves: u64,
    pub steps: u64,
    pub accepted: u64,
    pub rejected: u64,
    pub nonlinear_iterations: u64,
    pub succeeded: bool,
}

impl CostPoint {
    fn from(rtol: f64, atol: f64, error: f64, stats: &Stats, succeeded: bool) -> CostPoint {
        CostPoint {
            rtol,
            atol,
            error,
            rhs_evals: stats.rhs_evals,
            jacobian_evals: stats.jacobian_evals,
            lu_decompositions: stats.lu_decompositions,
            linear_solves: stats.linear_solves,
            steps: stats.steps,
            accepted: stats.accepted,
            rejected: stats.rejected,
            nonlinear_iterations: stats.nonlinear_iterations,
            succeeded,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct WorkPrecision {
    pub method: String,
    pub problem: String,
    pub points: Vec<CostPoint>,
    /// True when the errors were measured against a closed form solution.
    pub exact_reference: bool,
}

fn relative_error(computed: &[f64], reference: &[f64]) -> f64 {
    let mut worst = 0.0f64;
    for i in 0..reference.len() {
        let scale = reference[i].abs().max(1e-12);
        worst = worst.max((computed[i] - reference[i]).abs() / scale);
    }
    worst
}

/// The usual tolerance ladder, decade by decade.
pub fn tolerance_ladder(from: i32, to: i32) -> Vec<f64> {
    (to..=from).rev().map(|e| 10f64.powi(e)).collect()
}

/// Run one method over a tolerance ladder on one problem.
pub fn work_precision(
    method: &Method,
    problem: &dyn TestProblem,
    tolerances: &[f64],
    reference_method: Option<&Method>,
    template: &Options,
) -> WorkPrecision {
    let t_span = problem.t_span();
    let y0 = problem.y0();

    let (target, exact_reference) = match problem.exact(t_span.1) {
        Some(exact) => (Some(exact), true),
        None => {
            let reference = reference_method
                .and_then(|m| super::convergence::reference_solution(m, problem, t_span, &y0));
            (reference, false)
        }
    };

    let mut points = Vec::new();
    if let Some(target) = &target {
        for &rtol in tolerances {
            // The absolute tolerance follows the relative one by the usual
            // three decade offset, which keeps small components meaningful.
            let atol = rtol * 1e-3;
            let mut options = template.clone();
            options.rtol = rtol;
            options.atol = atol;
            options.adaptive = true;
            options.h0 = None;

            let solution = ode::integrate(method, problem, t_span, &y0, &options);
            let error = match solution.last() {
                Some(y) if solution.succeeded() => relative_error(y, target),
                _ => f64::INFINITY,
            };
            points.push(CostPoint::from(
                rtol,
                atol,
                error,
                &solution.stats,
                solution.succeeded(),
            ));
        }
    }

    WorkPrecision {
        method: method.id.clone(),
        problem: problem.id().to_string(),
        points,
        exact_reference,
    }
}

/// A single scalar cost, for ranking methods.
///
/// `jacobian_weight` and `factorization_weight` say how expensive a Jacobian
/// and a factorization are relative to one right hand side evaluation, which
/// is what turns the component counts into a comparable number.
pub fn weighted_cost(
    point: &CostPoint,
    dim: usize,
    jacobian_weight: f64,
    factorization_weight: f64,
) -> f64 {
    let n = dim as f64;
    point.rhs_evals as f64
        + jacobian_weight * n * point.jacobian_evals as f64
        + factorization_weight * n * n * n / 3.0 * point.lu_decompositions as f64
        + n * n * point.linear_solves as f64
}

/// Reference implementation of the work metric a stiff solver is judged by.
pub fn default_cost(point: &CostPoint, dim: usize) -> f64 {
    weighted_cost(point, dim, 1.0, 1.0)
}

/// Fraction of attempted steps that were thrown away, a direct read on how well
/// the controller is tuned for the problem.
pub fn rejection_rate(point: &CostPoint) -> f64 {
    if point.steps == 0 {
        return 0.0;
    }
    point.rejected as f64 / point.steps as f64
}

/// Convenience wrapper used by the command line tool and the browser build.
pub fn run<P: Problem + ?Sized>(
    method: &Method,
    problem: &P,
    t_span: (f64, f64),
    y0: &[f64],
    options: &Options,
) -> (Vec<f64>, Stats, bool) {
    let solution = ode::integrate(method, problem, t_span, y0, options);
    let succeeded = solution.succeeded();
    let last = solution.last().cloned().unwrap_or_default();
    (last, solution.stats, succeeded)
}
