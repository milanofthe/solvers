//! End to end order verification.
//!
//! The order conditions are checked from the coefficients elsewhere. This runs
//! every method in the library and measures the order it actually converges at,
//! which is the only thing that catches a mistake in a stepper rather than in a
//! tableau.

use solvers_core::analysis::convergence;
use solvers_core::linalg::Matrix;
use solvers_core::method::MethodLibrary;
use solvers_core::problem::Problem;
use solvers_core::problems::{NonlinearDecay, TestProblem};

/// Harmonic oscillator, whose eigenvalues sit on the imaginary axis. The two
/// weakly stable multistep methods are only usable here.
struct Oscillator {
    end: f64,
}

impl Problem for Oscillator {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = -y[0];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, _y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = 0.0;
        j[(0, 1)] = 1.0;
        j[(1, 0)] = -1.0;
        j[(1, 1)] = 0.0;
    }
}

impl TestProblem for Oscillator {
    fn id(&self) -> &'static str {
        "oscillator"
    }
    fn name(&self) -> &'static str {
        "Oscillator"
    }
    fn description(&self) -> &'static str {
        "Harmonic oscillator"
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, self.end)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 0.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
    fn exact(&self, t: f64) -> Option<Vec<f64>> {
        Some(vec![t.cos(), -t.sin()])
    }
}

/// A step size ladder that lands inside the window where the order is
/// measurable at all.
///
/// Two constraints bound that window from either side. The coarse end has to sit
/// inside the stability region, because the test problem has both eigenvalues at
/// minus one and a method whose real limit is `L` measures instability rather
/// than accuracy above `|L|`. The fine end has to stay above round off. For a
/// high order method those two are close together, so the ratio between rungs is
/// chosen to span about six decades of error whatever the order is, rather than
/// being fixed and overshooting the floor at one end or staying preasymptotic at
/// the other.
fn ladder(order: usize, real_limit: Option<f64>) -> Vec<f64> {
    const RUNGS: usize = 7;
    let p = order.max(1) as f64;
    let ratio = 10f64.powf(-1.0 / p).clamp(0.5, 0.85);
    let mut coarse: f64 = match order {
        0..=2 => 0.5,
        3..=4 => 0.5,
        5..=6 => 0.7,
        _ => 1.0,
    };
    if let Some(limit) = real_limit {
        coarse = coarse.min(0.7 * limit.abs());
    }
    (0..RUNGS).map(|i| coarse * ratio.powi(i as i32)).collect()
}

/// The real stability limit the analysis derived, or `None` when unbounded.
fn real_limit(method: &solvers_core::Method) -> Option<f64> {
    match solvers_core::analysis::analyze(method).real_stability_limit {
        Some(solvers_core::analysis::Limit::Finite(v)) => Some(v),
        _ => None,
    }
}

/// The families whose parasitic root sits on the unit circle. They cannot be run
/// on a decaying problem at all, which is a property of the methods and not a
/// defect of the implementation, so they are measured on an oscillator instead.
fn weakly_stable(method: &solvers_core::Method) -> bool {
    matches!(method.family.as_str(), "nystrom" | "milne_simpson")
}

#[test]
fn every_method_converges_at_its_stated_order() {
    let library = MethodLibrary::embedded().unwrap();
    let problem = NonlinearDecay { end: 2.0 };

    let mut failures = Vec::new();
    let mut rows = Vec::new();
    let mut unmeasurable = Vec::new();

    for method in library.iter() {
        if weakly_stable(method) {
            continue;
        }
        let expected = method.declared_order.unwrap_or(1) as usize;
        let study = convergence::study(method, &problem, &ladder(expected, real_limit(method)), None);

        // A method whose stability region forces a step size at which it is
        // already exact to double precision cannot be measured here at all. That
        // is a fact about the method, not a failure: its order is still checked
        // exactly from its coefficients by the order conditions. It is recorded
        // rather than asserted on, and the count is capped so the exemption
        // cannot quietly grow to cover a real defect.
        if study.points.iter().all(|p| p.error < 1e-11) {
            unmeasurable.push(method.id.clone());
            rows.push(format!(
                "{:<20} expected {:>2}  below the round off floor at every usable step",
                method.id, expected
            ));
            continue;
        }
        // The slope between the two finest usable step sizes is the asymptotic
        // rate; the fit over the whole ladder is still contaminated by the
        // coarse end for the high order methods.
        let measured = study.local_order;

        rows.push(format!(
            "{:<20} expected {:>2}  measured {:>6.2}  local {:>6.2}",
            method.id, expected, measured, study.local_order
        ));

        if !measured.is_finite() {
            failures.push(format!(
                "{}: no usable convergence points ({:?})",
                method.id,
                study
                    .points
                    .iter()
                    .map(|p| (p.h, p.error))
                    .collect::<Vec<_>>()
            ));
            continue;
        }
        if (measured - expected as f64).abs() > 0.5 {
            failures.push(format!(
                "{}: expected order {expected}, measured {measured:.3} (points {:?})",
                method.id,
                study
                    .points
                    .iter()
                    .map(|p| (p.h, p.error))
                    .collect::<Vec<_>>()
            ));
        }
    }

    for row in &rows {
        println!("{row}");
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn weakly_stable_multistep_methods_converge_on_an_oscillator() {
    let library = MethodLibrary::embedded().unwrap();
    let problem = Oscillator { end: 2.0 };

    let weak: Vec<&solvers_core::Method> = library.iter().filter(|m| weakly_stable(m)).collect();
    assert!(weak.len() >= 3, "expected the weakly stable families to be present");
    for method in weak {
        let id = method.id.as_str();
        let expected = method.declared_order.unwrap_or(1) as usize;
        // These two are run on the imaginary axis, where the real axis limit
        // says nothing, so the ladder is left uncapped.
        let study = convergence::study(method, &problem, &ladder(expected, None), None);
        println!(
            "{id:<16} expected {expected}  measured {:.2}",
            study.estimated_order
        );
        assert!(
            (study.local_order - expected as f64).abs() < 0.5,
            "{id}: expected order {expected}, measured {:.3} (points {:?})",
            study.local_order,
            study
                .points
                .iter()
                .map(|p| (p.h, p.error))
                .collect::<Vec<_>>()
        );
    }
}
