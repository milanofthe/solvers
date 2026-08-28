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
use solvers_core::problems::TestProblem;

/// `y1' = -y1 + y2^2`, `y2' = -y2`, with the closed form solution
/// `(2 exp(-t) - exp(-2t), exp(-t))`.
///
/// A linear problem would be the obvious order probe, but it only measures the
/// order the stability function agrees with the exponential to, and for several
/// high order methods that is higher than their true order. The quadratic term
/// here makes the elementary differentials nonzero so the full set of tree
/// conditions is exercised, while both eigenvalues stay at minus one so even
/// the methods with a small stability region have room to work.
struct Decay {
    end: f64,
}

impl Problem for Decay {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -y[0] + y[1] * y[1];
        dy[1] = -y[1];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = -1.0;
        j[(0, 1)] = 2.0 * y[1];
        j[(1, 0)] = 0.0;
        j[(1, 1)] = -1.0;
    }
}

impl TestProblem for Decay {
    fn id(&self) -> &'static str {
        "decay"
    }
    fn name(&self) -> &'static str {
        "Decay"
    }
    fn description(&self) -> &'static str {
        "Scalar decay"
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, self.end)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 1.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
    fn exact(&self, t: f64) -> Option<Vec<f64>> {
        Some(vec![2.0 * (-t).exp() - (-2.0 * t).exp(), (-t).exp()])
    }
}

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

/// Step sizes chosen so the errors land between the preasymptotic range and
/// round off, which for a high order method means a much coarser ladder.
///
/// The coarse end is also held inside the stability region: the test problem
/// has both eigenvalues at minus one, so a method whose real stability limit is
/// `L` must be run with `h` well below `|L|` or the coarse points measure
/// instability instead of accuracy.
fn ladder(order: usize, real_limit: Option<f64>) -> Vec<f64> {
    let (mut coarse, ratio): (f64, f64) = match order {
        0..=2 => (0.5, 0.5),
        3..=4 => (0.5, 0.6),
        5..=6 => (0.7, 0.65),
        _ => (1.0, 0.7),
    };
    if let Some(limit) = real_limit {
        coarse = coarse.min(0.4 * limit.abs());
    }
    (0..7).map(|i| coarse * ratio.powi(i)).collect()
}

/// The real stability limit the analysis derived, or `None` when unbounded.
fn real_limit(method: &solvers_core::Method) -> Option<f64> {
    match solvers_core::analysis::analyze(method).real_stability_limit {
        Some(solvers_core::analysis::Limit::Finite(v)) => Some(v),
        _ => None,
    }
}

/// The two methods whose parasitic root sits on the unit circle. They cannot be
/// run on a decaying problem at all, which is a property of the methods and not
/// a defect of the implementation.
const WEAKLY_STABLE: [&str; 2] = ["nystrom2", "milne_simpson"];

#[test]
fn every_method_converges_at_its_stated_order() {
    let library = MethodLibrary::embedded().unwrap();
    let problem = Decay { end: 2.0 };

    let mut failures = Vec::new();
    let mut rows = Vec::new();

    for method in library.iter() {
        if WEAKLY_STABLE.contains(&method.id.as_str()) {
            continue;
        }
        let expected = method.declared_order.unwrap_or(1) as usize;
        let study = convergence::study(method, &problem, &ladder(expected, real_limit(method)), None);
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

    for id in WEAKLY_STABLE {
        let method = library.get(id).unwrap();
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
