//! Command line access to the method library and its analyses.
//!
//! The same entry points the browser build exposes, for scripting and for
//! checking the library in continuous integration.

use solvers_core::analysis::{self, convergence, cost};
use solvers_core::method::{Method, MethodLibrary};
use solvers_core::ode::{self, Options};
use solvers_core::problems::{self, TestProblem};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let command = args.first().map(String::as_str).unwrap_or("help");

    let library = match MethodLibrary::embedded() {
        Ok(library) => library,
        Err(error) => {
            eprintln!("method library failed to load: {error}");
            std::process::exit(1);
        }
    };

    match command {
        "list" => list(&library),
        "problems" => list_problems(),
        "analyze" => analyze(&library, args.get(1)),
        "verify" => verify(&library),
        "convergence" => convergence_command(&library, &args),
        "cost" => cost_command(&library, &args),
        "solve" => solve_command(&library, &args),
        _ => usage(),
    }
}

fn usage() {
    println!("usage: solvers <command> [arguments]");
    println!();
    println!("  list                              every method in the library");
    println!("  problems                          every test problem");
    println!("  analyze <method>                  order and stability from the coefficients, as JSON");
    println!("  verify                            check every method file against its analysis");
    println!("  convergence <method> <problem>    fixed step convergence study, as JSON");
    println!("  cost <method> <problem>           adaptive work precision diagram, as JSON");
    println!("  solve <method> <problem> [rtol]   integrate and report the statistics");
}

fn find<'a>(library: &'a MethodLibrary, id: Option<&String>) -> &'a Method {
    let Some(id) = id else {
        eprintln!("a method id is required, see `solvers list`");
        std::process::exit(2);
    };
    match library.get(id) {
        Some(method) => method,
        None => {
            eprintln!("unknown method: {id}");
            std::process::exit(2);
        }
    }
}

fn find_problem(id: Option<&String>) -> Box<dyn TestProblem> {
    let Some(id) = id else {
        eprintln!("a problem id is required, see `solvers problems`");
        std::process::exit(2);
    };
    match problems::get(id) {
        Some(problem) => problem,
        None => {
            eprintln!("unknown problem: {id}");
            std::process::exit(2);
        }
    }
}

fn list(library: &MethodLibrary) {
    let mut rows: Vec<&Method> = library.iter().collect();
    rows.sort_by(|a, b| (a.family.as_str(), a.id.as_str()).cmp(&(&b.family, &b.id)));
    println!(
        "{:<22} {:<16} {:>5} {:>7} {:>10}  {}",
        "id", "family", "order", "size", "stability", "name"
    );
    for method in rows {
        let report = analysis::analyze(method);
        let stability = if report.l_stable {
            "L-stable".to_string()
        } else if report.a_stable {
            "A-stable".to_string()
        } else if let Some(angle) = report.alpha_angle {
            format!("A({angle:.1})")
        } else {
            "conditional".to_string()
        };
        println!(
            "{:<22} {:<16} {:>5} {:>7} {:>10}  {}",
            method.id, method.family, report.computed_order, report.size, stability, method.name
        );
    }
    println!("\n{} methods", library.len());
}

fn list_problems() {
    println!("{:<20} {:>4} {:>7}  {}", "id", "dim", "stiff", "name");
    for problem in problems::catalog() {
        println!(
            "{:<20} {:>4} {:>7}  {}",
            problem.id(),
            problem.dim(),
            if problem.is_stiff() { "yes" } else { "no" },
            problem.name()
        );
    }
}

fn analyze(library: &MethodLibrary, id: Option<&String>) {
    let method = find(library, id);
    let report = analysis::analyze(method);
    println!("{}", serde_json::to_string_pretty(&report).unwrap());
}

fn verify(library: &MethodLibrary) {
    let mut failures = 0;
    for method in library.iter() {
        let report = analysis::analyze(method);
        if report.discrepancies.is_empty() {
            continue;
        }
        failures += 1;
        println!("{}:", method.id);
        for issue in &report.discrepancies {
            println!("  {issue}");
        }
    }
    if failures == 0 {
        println!("all {} methods agree with their analysis", library.len());
    } else {
        println!("\n{failures} of {} methods disagree", library.len());
        std::process::exit(1);
    }
}

/// A high accuracy method for reference solutions.
fn reference<'a>(library: &'a MethodLibrary, stiff: bool) -> Option<&'a Method> {
    library.get(if stiff { "esdirk85" } else { "rkdp87" })
}

fn convergence_command(library: &MethodLibrary, args: &[String]) {
    let method = find(library, args.get(1));
    let problem = find_problem(args.get(2));
    let order = analysis::analyze(method).computed_order;
    let coarse = if order >= 7 { 0.5 } else if order >= 5 { 0.4 } else { 0.3 };
    let ratio: f64 = if order >= 5 { 0.7 } else { 0.6 };
    let steps: Vec<f64> = (0..7).map(|i| coarse * ratio.powi(i)).collect();
    let study = convergence::study(
        method,
        problem.as_ref(),
        &steps,
        reference(library, problem.is_stiff()),
    );
    println!("{}", serde_json::to_string_pretty(&study).unwrap());
}

fn cost_command(library: &MethodLibrary, args: &[String]) {
    let method = find(library, args.get(1));
    let problem = find_problem(args.get(2));
    let tolerances = cost::tolerance_ladder(-3, -11);
    let mut template = Options::default();
    template.max_steps = 2_000_000;
    let result = cost::work_precision(
        method,
        problem.as_ref(),
        &tolerances,
        reference(library, problem.is_stiff()),
        &template,
    );
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

fn solve_command(library: &MethodLibrary, args: &[String]) {
    let method = find(library, args.get(1));
    let problem = find_problem(args.get(2));
    let rtol: f64 = args
        .get(3)
        .and_then(|v| v.parse().ok())
        .unwrap_or(1e-6);

    let mut options = Options::with_tolerances(rtol, rtol * 1e-2);
    options.max_steps = 5_000_000;
    let span = problem.t_span();
    let y0 = problem.y0();

    let start = std::time::Instant::now();
    let solution = ode::integrate(method, problem.as_ref(), span, &y0, &options);
    let elapsed = start.elapsed();

    println!("{} on {}", method.name, problem.name());
    println!("  interval        {} to {}", span.0, span.1);
    println!("  status          {:?}", solution.status);
    println!("  wall time       {elapsed:.1?}");
    println!("  accepted steps  {}", solution.stats.accepted);
    println!("  rejected steps  {}", solution.stats.rejected);
    println!("  rhs evaluations {}", solution.stats.rhs_evals);
    println!("  jacobians       {}", solution.stats.jacobian_evals);
    println!("  factorizations  {}", solution.stats.lu_decompositions);
    println!("  newton steps    {}", solution.stats.nonlinear_iterations);

    if let (Some(y), Some(exact)) = (solution.last(), problem.exact(span.1)) {
        let error = y
            .iter()
            .zip(&exact)
            .fold(0.0f64, |acc, (a, b)| acc.max((a - b).abs() / b.abs().max(1e-10)));
        println!("  relative error  {error:.3e}");
    } else if let Some(y) = solution.last() {
        let formatted: Vec<String> = y.iter().map(|v| format!("{v:.6e}")).collect();
        println!("  final state     [{}]", formatted.join(", "));
    }
}
