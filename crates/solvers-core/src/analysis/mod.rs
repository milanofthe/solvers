//! Automated analysis of a method from its coefficients alone.
//!
//! Nothing here reads the claimed properties of a method file except to compare
//! against them. Order, stage order, stability and damping are all derived, so
//! a typo in a tableau shows up as a disagreement rather than as a method that
//! quietly integrates at the wrong order.

pub mod convergence;
pub mod cost;
pub mod order;
pub mod stability;

pub use order::{OrderReport, Tree};
pub use stability::{GeneratingPolynomials, StabilityFunction, StabilityGrid};

use crate::method::{Method, MethodKind};
use crate::num::Complex;
use serde::Serialize;

/// A stability limit along a ray, which may legitimately be infinite.
///
/// JSON has no infinity, and reporting a huge number where the truth is "the
/// whole ray is stable" would be read as a real bound, so the unbounded case
/// gets its own value.
#[derive(Copy, Clone, Debug, Serialize)]
#[serde(untagged)]
pub enum Limit {
    Finite(f64),
    Unbounded(Unbounded),
}

#[derive(Copy, Clone, Debug, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Unbounded {
    Unbounded,
}

impl From<f64> for Limit {
    fn from(value: f64) -> Limit {
        if value.is_finite() {
            Limit::Finite(value)
        } else {
            Limit::Unbounded(Unbounded::Unbounded)
        }
    }
}

/// Everything that can be said about a method without running it.
#[derive(Clone, Debug, Serialize)]
pub struct MethodReport {
    pub id: String,
    pub name: String,
    pub family: String,
    pub class: &'static str,
    /// Stages for Runge-Kutta, steps for multistep.
    pub size: usize,
    pub implicit: bool,
    pub adaptive: bool,

    pub declared_order: Option<u32>,
    pub computed_order: usize,
    pub declared_embedded_order: Option<u32>,
    pub computed_embedded_order: Option<usize>,
    pub stage_order: Option<usize>,
    pub consistent_abscissae: Option<bool>,
    /// True when the order conditions were checked in exact arithmetic.
    pub exact_arithmetic: bool,

    pub a_stable: bool,
    pub l_stable: bool,
    pub stiffly_accurate: Option<bool>,
    /// `R(infinity)` for Runge-Kutta methods.
    pub damping_at_infinity: Option<f64>,
    /// Half angle of the A(alpha) wedge in degrees, ninety meaning A-stable.
    pub alpha_angle: Option<f64>,
    pub real_stability_limit: Option<Limit>,
    pub imaginary_stability_limit: Option<Limit>,
    pub zero_stable: Option<bool>,

    /// Effective cost of one step in right hand side evaluations.
    pub stage_cost: usize,
    /// Points where the method file and the analysis disagree.
    pub discrepancies: Vec<String>,
}

/// Derive everything about a method from its coefficients.
pub fn analyze(method: &Method) -> MethodReport {
    let mut discrepancies = Vec::new();

    match &method.kind {
        MethodKind::RungeKutta(tableau) => {
            let report = order::verify(tableau, 10);
            let function = StabilityFunction::from_tableau(tableau);
            let a_stable = function.is_a_stable();
            let l_stable = function.is_l_stable();

            if let Some(declared) = method.declared_order {
                if declared as usize != report.order {
                    discrepancies.push(format!(
                        "file claims order {declared}, the tableau satisfies order {}",
                        report.order
                    ));
                }
            }
            if let (Some(declared), Some(computed)) =
                (method.declared_embedded_order, report.embedded_order)
            {
                if declared as usize != computed {
                    discrepancies.push(format!(
                        "file claims embedded order {declared}, the tableau satisfies {computed}"
                    ));
                }
            }
            if let Some(claimed) = method.properties.a_stable {
                if claimed != a_stable {
                    discrepancies.push(format!("file claims a_stable = {claimed}, computed {a_stable}"));
                }
            }
            if let Some(claimed) = method.properties.l_stable {
                if claimed != l_stable {
                    discrepancies.push(format!("file claims l_stable = {claimed}, computed {l_stable}"));
                }
            }
            if let Some(claimed) = method.properties.stiffly_accurate {
                if claimed != tableau.stiffly_accurate {
                    discrepancies.push(format!(
                        "file claims stiffly_accurate = {claimed}, computed {}",
                        tableau.stiffly_accurate
                    ));
                }
            }
            if let Some(claimed) = method.properties.stage_order {
                if claimed as usize != report.stage_order {
                    discrepancies.push(format!(
                        "file claims stage order {claimed}, computed {}",
                        report.stage_order
                    ));
                }
            }
            if !report.consistent_abscissae {
                discrepancies.push("c is not the row sum of A".to_string());
            }

            // Cost of a step: implicit stages need a solve, explicit ones one
            // evaluation. FSAL saves one evaluation per step.
            let explicit_stages = (0..tableau.stages)
                .filter(|&i| tableau.a[(i, i)].value().abs() < 1e-14)
                .count();
            let stage_cost = if tableau.is_fsal() {
                tableau.stages - 1
            } else {
                tableau.stages
            }
            .max(explicit_stages);

            MethodReport {
                id: method.id.clone(),
                name: method.name.clone(),
                family: method.family.clone(),
                class: "runge_kutta",
                size: tableau.stages,
                implicit: !tableau.is_explicit(),
                adaptive: tableau.has_embedded(),
                declared_order: method.declared_order,
                computed_order: report.order,
                declared_embedded_order: method.declared_embedded_order,
                computed_embedded_order: report.embedded_order,
                stage_order: Some(report.stage_order),
                consistent_abscissae: Some(report.consistent_abscissae),
                exact_arithmetic: report.exact,
                a_stable,
                l_stable,
                stiffly_accurate: Some(tableau.stiffly_accurate),
                damping_at_infinity: Some(function.at_infinity()),
                alpha_angle: if a_stable { Some(90.0) } else { None },
                real_stability_limit: Some(function.real_stability_limit().into()),
                imaginary_stability_limit: Some(function.imaginary_stability_limit().into()),
                zero_stable: None,
                stage_cost,
                discrepancies,
            }
        }
        MethodKind::LinearMultistep(family) => {
            let coefficients = family.uniform_coefficients();
            let (computed_order, polynomials) = match &coefficients {
                Ok(c) => (c.order, Some(GeneratingPolynomials::from_coefficients(c))),
                Err(e) => {
                    discrepancies.push(format!("coefficients could not be determined: {e}"));
                    (0, None)
                }
            };

            if let Some(declared) = method.declared_order {
                if declared as usize != computed_order {
                    discrepancies.push(format!(
                        "file claims order {declared}, the coefficients give order {computed_order}"
                    ));
                }
            }

            let alpha_angle = polynomials.as_ref().map(|p| p.alpha_angle());
            let zero_stable = polynomials.as_ref().map(|p| p.is_zero_stable());
            let a_stable = alpha_angle.map_or(false, |a| a >= 89.5);

            if let Some(claimed) = method.properties.a_stable {
                if claimed != a_stable {
                    discrepancies.push(format!("file claims a_stable = {claimed}, computed {a_stable}"));
                }
            }

            let real_limit = polynomials
                .as_ref()
                .map(|p| Limit::from(stability::scan_real_limit(|x| p.is_stable_at(Complex::real(x)))));

            MethodReport {
                id: method.id.clone(),
                name: method.name.clone(),
                family: method.family.clone(),
                class: "linear_multistep",
                size: family.steps,
                implicit: family.implicit,
                adaptive: family.steps > 1,
                declared_order: method.declared_order,
                computed_order,
                declared_embedded_order: method.declared_embedded_order,
                computed_embedded_order: None,
                stage_order: None,
                consistent_abscissae: None,
                exact_arithmetic: false,
                a_stable,
                l_stable: false,
                stiffly_accurate: None,
                damping_at_infinity: None,
                alpha_angle,
                real_stability_limit: real_limit,
                imaginary_stability_limit: None,
                zero_stable,
                // One implicit solve, or one evaluation for an explicit family.
                stage_cost: 1,
                discrepancies,
            }
        }
    }
}

/// Sample the stability region of any method on a rectangle.
pub fn stability_region(
    method: &Method,
    re: (f64, f64),
    im: (f64, f64),
    width: usize,
    height: usize,
) -> Option<StabilityGrid> {
    match &method.kind {
        MethodKind::RungeKutta(tableau) => {
            let function = StabilityFunction::from_tableau(tableau);
            Some(stability::sample_region(
                |z| function.eval(z).abs(),
                re,
                im,
                width,
                height,
            ))
        }
        MethodKind::LinearMultistep(family) => {
            let coefficients = family.uniform_coefficients().ok()?;
            let polynomials = GeneratingPolynomials::from_coefficients(&coefficients);
            // For a multistep method the natural indicator is the largest root
            // modulus, which plays the role `|R(z)|` plays for Runge-Kutta.
            Some(stability::sample_region(
                |z| polynomials.root_radius(z),
                re,
                im,
                width,
                height,
            ))
        }
    }
}
