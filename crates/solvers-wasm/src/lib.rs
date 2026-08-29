//! Browser bindings.
//!
//! The web interface runs the same code the command line does. Nothing is
//! precomputed and shipped as data: stability regions, order reports and work
//! precision diagrams are all evaluated in the page from the same method files,
//! so a change to a tableau cannot leave the plots showing something else.

use solvers_core::analysis::{self, convergence, cost, stability};
use solvers_core::control::ControllerPreset;
use solvers_core::method::{Method, MethodClass, MethodKind, MethodLibrary};
use solvers_core::nonlinear::SolverKind;
use solvers_core::num::Complex;
use solvers_core::ode::{self, Options};
use solvers_core::problems::{self, TestProblem};
use serde_json::{json, Value};
use std::sync::OnceLock;
use wasm_bindgen::prelude::*;

fn library() -> &'static MethodLibrary {
    static LIBRARY: OnceLock<MethodLibrary> = OnceLock::new();
    LIBRARY.get_or_init(|| MethodLibrary::embedded().expect("embedded method library must be valid"))
}

fn find(id: &str) -> Result<&'static Method, JsValue> {
    library()
        .get(id)
        .ok_or_else(|| JsValue::from_str(&format!("unknown method: {id}")))
}

fn find_problem(id: &str) -> Result<Box<dyn TestProblem>, JsValue> {
    problems::get(id).ok_or_else(|| JsValue::from_str(&format!("unknown problem: {id}")))
}

/// A high accuracy method to generate reference solutions with.
fn reference_method(stiff: bool) -> &'static Method {
    let id = if stiff { "esdirk85" } else { "rkdp87" };
    library().get(id).expect("reference method must exist")
}

/// Every method with the properties the catalogue view needs.
#[wasm_bindgen]
pub fn method_catalog() -> String {
    let entries: Vec<Value> = library()
        .iter()
        .map(|method| {
            let report = analysis::analyze(method);
            json!({
                "id": method.id,
                "name": method.name,
                "family": method.family,
                "class": report.class,
                "description": method.description,
                "size": report.size,
                "implicit": report.implicit,
                "adaptive": report.adaptive,
                "order": report.computed_order,
                "embeddedOrder": report.computed_embedded_order,
                "stageOrder": report.stage_order,
                "aStable": report.a_stable,
                "lStable": report.l_stable,
                "stifflyAccurate": report.stiffly_accurate,
                "alphaAngle": report.alpha_angle,
                "dampingAtInfinity": report.damping_at_infinity,
                "stageCost": report.stage_cost,
                "exactArithmetic": report.exact_arithmetic,
                "discrepancies": report.discrepancies,
                "doi": method.references.first().and_then(|r| r.doi.clone()),
            })
        })
        .collect();
    serde_json::to_string(&entries).unwrap_or_else(|_| "[]".into())
}

/// Full detail for one method: coefficients, analysis and references.
#[wasm_bindgen]
pub fn method_detail(id: &str) -> Result<String, JsValue> {
    let method = find(id)?;
    let report = analysis::analyze(method);

    let coefficients = match &method.kind {
        MethodKind::RungeKutta(tableau) => {
            let row = |values: &[solvers_core::num::Coeff]| -> Vec<Value> {
                values
                    .iter()
                    .map(|c| json!({ "text": c.to_string(), "value": c.value(), "exact": c.is_exact() }))
                    .collect()
            };
            let a: Vec<Vec<Value>> = (0..tableau.stages)
                .map(|i| row(&(0..tableau.stages).map(|j| tableau.a[(i, j)]).collect::<Vec<_>>()))
                .collect();
            json!({
                "kind": "runge_kutta",
                "stages": tableau.stages,
                "a": a,
                "b": row(&tableau.b),
                "c": row(&tableau.c),
                "bEmbedded": tableau.b_embedded.as_ref().map(|v| row(v)),
                "structure": tableau.structure,
                "singlyDiagonal": tableau.singly_diagonal,
                "explicitFirstStage": tableau.explicit_first_stage,
                "fsal": tableau.is_fsal(),
                "gamma": tableau.gamma.map(|g| g.to_string()),
            })
        }
        MethodKind::LinearMultistep(family) => {
            let coefficients = family.uniform_coefficients().ok();
            json!({
                "kind": "linear_multistep",
                "steps": family.steps,
                "minSteps": family.min_steps,
                "startup": family.startup,
                "alpha": coefficients.as_ref().map(|c| c.alpha.clone()),
                "beta": coefficients.as_ref().map(|c| c.beta.clone()),
            })
        }
    };

    let references: Vec<Value> = method
        .references
        .iter()
        .map(|r| {
            json!({
                "authors": r.authors,
                "title": r.title,
                "year": r.year,
                "source": r.source,
                "doi": r.doi,
                "link": r.link(),
            })
        })
        .collect();

    let detail = json!({
        "id": method.id,
        "name": method.name,
        "aliases": method.aliases,
        "family": method.family,
        "description": method.description,
        "report": report,
        "coefficients": coefficients,
        "references": references,
    });
    Ok(serde_json::to_string(&detail).unwrap_or_default())
}

/// `log10 |R(z)|` sampled on a rectangle, row major from `im_min` to `im_max`.
///
/// The log scale is what makes the plot readable: the interesting structure is
/// the unit contour, and around it the magnitude spans many decades.
#[wasm_bindgen]
pub fn stability_grid(
    id: &str,
    re_min: f64,
    re_max: f64,
    im_min: f64,
    im_max: f64,
    width: usize,
    height: usize,
) -> Result<Vec<f64>, JsValue> {
    let method = find(id)?;
    let grid = analysis::stability_region(method, (re_min, re_max), (im_min, im_max), width, height)
        .ok_or_else(|| JsValue::from_str("no stability region for this method"))?;
    Ok(grid
        .magnitude
        .into_iter()
        .map(|m| {
            if m.is_finite() && m > 0.0 {
                m.log10()
            } else if m == 0.0 {
                -30.0
            } else {
                30.0
            }
        })
        .collect())
}

/// The boundary locus of a multistep method, interleaved as `re, im, re, im`.
#[wasm_bindgen]
pub fn boundary_locus(id: &str, samples: usize) -> Result<Vec<f64>, JsValue> {
    let method = find(id)?;
    let family = method
        .multistep()
        .ok_or_else(|| JsValue::from_str("boundary locus is a multistep notion"))?;
    let coefficients = family
        .uniform_coefficients()
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let polynomials = stability::GeneratingPolynomials::from_coefficients(&coefficients);
    let mut out = Vec::with_capacity(2 * samples);
    for z in polynomials.boundary_locus(samples) {
        out.push(z.re);
        out.push(z.im);
    }
    Ok(out)
}

/// The stability function as a rational function, for display.
#[wasm_bindgen]
pub fn stability_function(id: &str) -> Result<String, JsValue> {
    let method = find(id)?;
    let tableau = method
        .tableau()
        .ok_or_else(|| JsValue::from_str("stability function needs a Runge-Kutta tableau"))?;
    let function = stability::StabilityFunction::from_tableau(tableau);
    let poles: Vec<Value> = function
        .poles()
        .iter()
        .map(|p| json!({ "re": p.re, "im": p.im }))
        .collect();
    let at_infinity = function.at_infinity();
    let value = json!({
        "numerator": function.numerator,
        "denominator": function.denominator,
        "atInfinity": if at_infinity.is_finite() {
            Value::from(at_infinity)
        } else {
            Value::from("unbounded")
        },
        "poles": poles,
        "orderOfConsistency": function.order_of_consistency(12),
        "atMinusOne": function.eval(Complex::real(-1.0)).abs(),
    });
    Ok(serde_json::to_string(&value).unwrap_or_default())
}

/// The test problems available to the analyses.
#[wasm_bindgen]
pub fn problem_catalog() -> String {
    let entries: Vec<Value> = problems::catalog()
        .iter()
        .map(|p| {
            let span = p.t_span();
            json!({
                "id": p.id(),
                "name": p.name(),
                "description": p.description(),
                "dim": p.dim(),
                "stiff": p.is_stiff(),
                "tStart": span.0,
                "tEnd": span.1,
                "hasExact": p.exact(span.0).is_some(),
            })
        })
        .collect();
    serde_json::to_string(&entries).unwrap_or_else(|_| "[]".into())
}

/// Fixed step convergence study.
#[wasm_bindgen]
pub fn convergence_study(
    method_id: &str,
    problem_id: &str,
    coarse: f64,
    ratio: f64,
    count: usize,
) -> Result<String, JsValue> {
    let method = find(method_id)?;
    let problem = find_problem(problem_id)?;
    let steps: Vec<f64> = (0..count.clamp(2, 20))
        .map(|i| coarse * ratio.powi(i as i32))
        .collect();
    let reference = reference_method(problem.is_stiff());
    let study = convergence::study(method, problem.as_ref(), &steps, Some(reference));
    Ok(serde_json::to_string(&study).unwrap_or_default())
}

/// Adaptive work precision diagram over a tolerance ladder.
#[wasm_bindgen]
pub fn work_precision(
    method_id: &str,
    problem_id: &str,
    from_exponent: i32,
    to_exponent: i32,
) -> Result<String, JsValue> {
    let method = find(method_id)?;
    let problem = find_problem(problem_id)?;
    let tolerances = cost::tolerance_ladder(from_exponent, to_exponent);
    let reference = reference_method(problem.is_stiff());
    let mut template = Options::default();
    template.max_steps = 2_000_000;
    let result = cost::work_precision(
        method,
        problem.as_ref(),
        &tolerances,
        Some(reference),
        &template,
    );
    Ok(serde_json::to_string(&result).unwrap_or_default())
}

/// One adaptive run, returned as a trajectory plus its statistics.
#[wasm_bindgen]
pub fn trajectory(
    method_id: &str,
    problem_id: &str,
    rtol: f64,
    atol: f64,
    samples: usize,
) -> Result<String, JsValue> {
    let method = find(method_id)?;
    let problem = find_problem(problem_id)?;
    let span = problem.t_span();
    let y0 = problem.y0();

    let mut options = Options::with_tolerances(rtol, atol);
    options.max_steps = 2_000_000;
    if samples > 1 {
        options.t_eval = Some(
            (0..samples)
                .map(|i| span.0 + (span.1 - span.0) * i as f64 / (samples - 1) as f64)
                .collect(),
        );
    }

    let solution = ode::integrate(method, problem.as_ref(), span, &y0, &options);
    let value = json!({
        "t": solution.t,
        "y": solution.y,
        "steps": solution.steps,
        "status": solution.status,
        "stats": solution.stats,
    });
    Ok(serde_json::to_string(&value).unwrap_or_default())
}

/// The step size history of an adaptive run, which is what a controller is
/// actually judged on.
#[wasm_bindgen]
pub fn step_history(
    method_id: &str,
    problem_id: &str,
    rtol: f64,
    atol: f64,
    controller: &str,
) -> Result<String, JsValue> {
    let method = find(method_id)?;
    let problem = find_problem(problem_id)?;
    let span = problem.t_span();
    let y0 = problem.y0();

    let mut options = Options::with_tolerances(rtol, atol);
    options.max_steps = 2_000_000;
    if let Some(config) = solvers_core::control::ControllerConfig::from_name(controller) {
        options.controller = config;
    }

    let solution = ode::integrate(method, problem.as_ref(), span, &y0, &options);
    let mut times = Vec::with_capacity(solution.steps.len());
    let mut t = span.0;
    for h in &solution.steps {
        times.push(t);
        t += h;
    }
    let value = json!({
        "t": times,
        "h": solution.steps,
        "stats": solution.stats,
        "status": solution.status,
    });
    Ok(serde_json::to_string(&value).unwrap_or_default())
}

/// Names of the available error controllers and nonlinear solvers.
#[wasm_bindgen]
pub fn options_catalog() -> String {
    let controllers: Vec<&str> = ControllerPreset::all().iter().map(|p| p.name()).collect();
    let solvers: Vec<&str> = SolverKind::all().iter().map(|k| k.name()).collect();
    serde_json::to_string(&json!({
        "controllers": controllers,
        "nonlinearSolvers": solvers,
        "classes": [MethodClass::RungeKutta, MethodClass::LinearMultistep],
    }))
    .unwrap_or_default()
}
