//! Adaptive integration: does the error control actually deliver the accuracy
//! it was asked for, and do the stiff methods survive a stiff problem.

use solvers_core::analysis::analyze;
use solvers_core::control::{ControllerConfig, ControllerPreset};
use solvers_core::method::MethodLibrary;
use solvers_core::nonlinear::SolverKind;
use solvers_core::ode::{self, Options};
use solvers_core::problems::{self, TestProblem};

fn relative_error(computed: &[f64], reference: &[f64]) -> f64 {
    let mut worst = 0.0f64;
    for i in 0..reference.len() {
        worst = worst.max((computed[i] - reference[i]).abs() / reference[i].abs().max(1e-10));
    }
    worst
}

fn run(
    method: &solvers_core::Method,
    problem: &dyn TestProblem,
    options: &Options,
) -> (f64, solvers_core::Stats, bool) {
    let span = problem.t_span();
    let y0 = problem.y0();
    let solution = ode::integrate(method, problem, span, &y0, options);
    let exact = problem.exact(span.1).expect("test needs a closed form solution");
    let error = match solution.last() {
        Some(y) if solution.succeeded() => relative_error(y, &exact),
        _ => f64::INFINITY,
    };
    (error, solution.stats, solution.succeeded())
}

#[test]
fn adaptive_runs_deliver_the_requested_accuracy() {
    let library = MethodLibrary::embedded().unwrap();
    let problem = problems::NonlinearDecay::default();

    let mut failures = Vec::new();
    for method in library.iter() {
        if !method.is_adaptive() {
            continue;
        }
        let mut options = Options::with_tolerances(1e-8, 1e-10);
        options.max_steps = 500_000;
        let (error, stats, ok) = run(method, &problem, &options);

        println!(
            "{:<20} error {:>10.3e}  steps {:>6} (rejected {:>4})  rhs {:>7}",
            method.id, error, stats.accepted, stats.rejected, stats.rhs_evals
        );

        if !ok {
            failures.push(format!("{}: integration did not succeed", method.id));
        } else if !(error < 1e-4) {
            // A tolerance of 1e-8 should land far inside 1e-4; anything worse
            // means the error estimate is not measuring what it should.
            failures.push(format!("{}: error {error:.3e} at rtol 1e-8", method.id));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn tightening_the_tolerance_improves_the_error() {
    let library = MethodLibrary::embedded().unwrap();
    let problem = problems::NonlinearDecay::default();

    for id in ["rkdp54", "rkbs32", "esdirk43", "bdf4", "adams_moulton_4"] {
        let method = library.get(id).unwrap();
        let mut loose = Options::with_tolerances(1e-4, 1e-6);
        loose.max_steps = 500_000;
        let mut tight = Options::with_tolerances(1e-9, 1e-11);
        tight.max_steps = 500_000;

        let (loose_error, loose_stats, _) = run(method, &problem, &loose);
        let (tight_error, tight_stats, _) = run(method, &problem, &tight);
        println!(
            "{id:<18} loose {loose_error:.3e} ({} steps)  tight {tight_error:.3e} ({} steps)",
            loose_stats.accepted, tight_stats.accepted
        );
        assert!(
            tight_error < loose_error,
            "{id}: tightening the tolerance did not reduce the error"
        );
        assert!(
            tight_stats.accepted > loose_stats.accepted,
            "{id}: tightening the tolerance did not cost more steps"
        );
    }
}

#[test]
fn stiff_methods_survive_a_stiff_problem() {
    let library = MethodLibrary::embedded().unwrap();
    // Prothero-Robinson with lambda = -1e4. The solution is cos(t) regardless
    // of the stiffness, so the accuracy of a stiff method is directly readable.
    let problem = problems::ProtheroRobinson { lambda: -1e4 };

    let mut tested = 0;
    let mut failures = Vec::new();
    for method in library.iter() {
        let report = analyze(method);
        if !report.a_stable || !method.is_adaptive() {
            continue;
        }
        let mut options = Options::with_tolerances(1e-6, 1e-8);
        options.max_steps = 200_000;
        let (error, stats, ok) = run(method, &problem, &options);
        tested += 1;
        println!(
            "{:<20} error {:>10.3e}  steps {:>6}  jac {:>5}  lu {:>5}",
            method.id, error, stats.accepted, stats.jacobian_evals, stats.lu_decompositions
        );
        if !ok {
            failures.push(format!("{}: integration did not succeed", method.id));
        } else if !(error < 1e-2) {
            failures.push(format!("{}: error {error:.3e} on Prothero-Robinson", method.id));
        }
    }
    assert!(tested >= 5, "expected the library to contain A-stable adaptive methods");
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn every_controller_preset_integrates_correctly() {
    let library = MethodLibrary::embedded().unwrap();
    let method = library.get("rkdp54").unwrap();
    let problem = problems::NonlinearDecay::default();

    for preset in ControllerPreset::all() {
        let mut options = Options::with_tolerances(1e-8, 1e-10);
        options.controller = ControllerConfig::preset(*preset);
        options.max_steps = 500_000;
        let (error, stats, ok) = run(method, &problem, &options);
        println!(
            "{:<14} error {:>10.3e}  steps {:>6}  rejected {:>4}",
            preset.name(),
            error,
            stats.accepted,
            stats.rejected
        );
        assert!(ok, "{}: integration failed", preset.name());
        assert!(error < 1e-4, "{}: error {error:.3e}", preset.name());
    }
}

#[test]
fn every_nonlinear_solver_drives_an_implicit_method() {
    let library = MethodLibrary::embedded().unwrap();
    let method = library.get("esdirk43").unwrap();
    let problem = problems::NonlinearDecay::default();

    for kind in SolverKind::all() {
        let mut options = Options::with_tolerances(1e-8, 1e-10);
        options.nonlinear.kind = *kind;
        options.max_steps = 500_000;
        let (error, stats, ok) = run(method, &problem, &options);
        println!(
            "{:<18} error {:>10.3e}  steps {:>6}  iterations {:>7}",
            kind.name(),
            error,
            stats.accepted,
            stats.nonlinear_iterations
        );
        assert!(ok, "{}: integration failed", kind.name());
        assert!(error < 1e-4, "{}: error {error:.3e}", kind.name());
    }
}

#[test]
fn dense_output_matches_the_solution_on_the_grid() {
    let library = MethodLibrary::embedded().unwrap();
    let problem = problems::NonlinearDecay::default();
    let span = problem.t_span();
    let y0 = problem.y0();

    for id in ["rkdp54", "rkbs32", "esdirk43", "bdf4"] {
        let method = library.get(id).unwrap();
        let grid: Vec<f64> = (0..=20).map(|i| span.0 + (span.1 - span.0) * i as f64 / 20.0).collect();
        let mut options = Options::with_tolerances(1e-10, 1e-12);
        options.t_eval = Some(grid.clone());
        options.max_steps = 500_000;

        let solution = ode::integrate(method, &problem, span, &y0, &options);
        assert!(solution.succeeded(), "{id}: integration failed");
        assert_eq!(solution.t.len(), grid.len(), "{id}: wrong number of output points");

        let mut worst = 0.0f64;
        for (t, y) in solution.t.iter().zip(&solution.y) {
            let exact = problem.exact(*t).unwrap();
            worst = worst.max(relative_error(y, &exact));
        }
        println!("{id:<12} worst interpolation error {worst:.3e}");
        // Interpolation is allowed to be less accurate than the step points,
        // but not by orders of magnitude.
        assert!(worst < 1e-3, "{id}: interpolation error {worst:.3e}");
    }
}
