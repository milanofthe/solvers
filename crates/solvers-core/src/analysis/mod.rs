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
pub mod tags;

pub use order::{OrderReport, Tree};
pub use stability::{GeneratingPolynomials, StabilityFunction, StabilityGrid};
pub use tags::{tags, Tag};

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
    /// `R(infinity)` for Runge-Kutta methods. Unbounded for an explicit method,
    /// whose stability function is a polynomial.
    pub damping_at_infinity: Option<Limit>,
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
                damping_at_infinity: Some(function.at_infinity().into()),
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

/// A window on the complex plane that actually contains the picture.
///
/// A fixed window per method class is always wrong for somebody: backward Euler
/// needs a couple of units, Radau needs a dozen, and an eight stage explicit
/// pair needs more still. So the window is measured rather than guessed.
///
/// For a multistep method the boundary is known in closed form and there is
/// nothing to probe for: the curve is the frame. For a Runge-Kutta method a
/// coarse probe is sampled instead, and whichever of the stable and the
/// unstable set turns out to be bounded is the one that sets the extent; the
/// other is the half plane the method lives in and has no extent to speak of.
pub fn suggested_window(method: &Method, aspect: f64) -> ((f64, f64), (f64, f64)) {
    let scale = probe_scale(method);
    let fallback = shape_to_aspect((-scale * 0.6, scale * 0.4), (-scale * 0.5, scale * 0.5), aspect);

    let chosen = match boundary_bounds(method).or_else(|| probe_bounds(method, scale)) {
        Some(bounds) => bounds,
        None => return fallback,
    };

    // The origin is the reference point of every one of these pictures, so it
    // stays in view even when the interesting set sits away from it.
    let mut x0 = chosen.x0.min(0.0);
    let mut x1 = chosen.x1.max(0.0);
    let mut y0 = chosen.y0.min(0.0);
    let mut y1 = chosen.y1.max(0.0);

    // One margin for both axes, taken from the larger side. A fixed minimum
    // would swallow the picture whole for a family like Adams-Bashforth 8,
    // whose region is a fraction of a unit across.
    let pad = (x1 - x0).max(y1 - y0).max(1e-9) * 0.18;
    x0 -= pad;
    x1 += pad;
    y0 -= pad;
    y1 += pad;

    shape_to_aspect((x0, x1), (y0, y1), aspect)
}

/// The extent of a multistep stability region, taken from its own boundary.
///
/// Exact where a sampled grid is only as good as its resolution, and the only
/// thing that works at all for a region with no interior: Nystrom 2 and
/// Milne-Simpson are stable on a segment of the imaginary axis and nowhere
/// else, and no grid resolves a set of zero area.
///
/// Two cases have no box to give and hand the question back. The boundary of an
/// A-stable family runs to infinity, which shows up as a tail orders of
/// magnitude past the middle of the curve. And Nystrom above two steps is
/// stable at the origin and nowhere else, which leaves no extent at all.
fn boundary_bounds(method: &Method) -> Option<Bounds> {
    let family = method.multistep()?;
    let coefficients = family.uniform_coefficients().ok()?;
    let polynomials = GeneratingPolynomials::from_coefficients(&coefficients);

    let points: Vec<Complex> = polynomials
        .region_boundary(720)
        .into_iter()
        .filter(|z| z.abs().is_finite())
        .collect();
    if points.is_empty() {
        return None;
    }
    let mut moduli: Vec<f64> = points.iter().map(|z| z.abs()).collect();
    moduli.sort_by(f64::total_cmp);
    if *moduli.last().unwrap() > 1e3 * moduli[moduli.len() / 2].max(1e-12) {
        return None;
    }

    let mut bounds = Bounds::empty();
    for z in &points {
        bounds.include(z.re, z.im, false);
    }
    let extent = (bounds.x1 - bounds.x0).max(bounds.y1 - bounds.y0);
    (extent > 1e-9).then_some(bounds)
}

/// The extent of whichever of the two sets a coarse probe finds bounded.
fn probe_bounds(method: &Method, scale: f64) -> Option<Bounds> {
    let probe = (-scale, scale);
    let samples = 72;
    let grid = stability_region(method, probe, probe, samples, samples)?;

    let mut stable = Bounds::empty();
    let mut unstable = Bounds::empty();
    for row in 0..grid.height {
        let y = probe.0 + (probe.1 - probe.0) * row as f64 / (grid.height - 1) as f64;
        for column in 0..grid.width {
            let x = probe.0 + (probe.1 - probe.0) * column as f64 / (grid.width - 1) as f64;
            let value = grid.magnitude[row * grid.width + column];
            let edge = row == 0 || column == 0 || row + 1 == grid.height || column + 1 == grid.width;
            if value.is_finite() && value <= 1.0 {
                stable.include(x, y, edge);
            } else {
                unstable.include(x, y, edge);
            }
        }
    }

    // The bounded set is the one that does not run off the probe. Where neither
    // is, the boundary is a line through the plane, as it is for a method whose
    // region is exactly a half plane, and no box means anything.
    if stable.bounded() {
        Some(stable)
    } else if unstable.bounded() {
        Some(unstable)
    } else {
        None
    }
}

#[derive(Clone, Copy)]
struct Bounds {
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    touches_edge: bool,
    any: bool,
}

impl Bounds {
    fn empty() -> Bounds {
        Bounds {
            x0: f64::INFINITY,
            x1: f64::NEG_INFINITY,
            y0: f64::INFINITY,
            y1: f64::NEG_INFINITY,
            touches_edge: false,
            any: false,
        }
    }

    fn include(&mut self, x: f64, y: f64, edge: bool) {
        self.x0 = self.x0.min(x);
        self.x1 = self.x1.max(x);
        self.y0 = self.y0.min(y);
        self.y1 = self.y1.max(y);
        self.touches_edge |= edge;
        self.any = true;
    }

    fn bounded(&self) -> bool {
        self.any && !self.touches_edge
    }
}

/// Grow the shorter side until the box has the aspect the panel does, so the
/// data fills the frame without a shape being stretched out of true.
fn shape_to_aspect(re: (f64, f64), im: (f64, f64), aspect: f64) -> ((f64, f64), (f64, f64)) {
    let width = re.1 - re.0;
    let height = im.1 - im.0;
    if width <= 0.0 || height <= 0.0 || !aspect.is_finite() || aspect <= 0.0 {
        return (re, im);
    }
    if width / height < aspect {
        let target = height * aspect;
        let centre = 0.5 * (re.0 + re.1);
        ((centre - target / 2.0, centre + target / 2.0), im)
    } else {
        let target = width / aspect;
        let centre = 0.5 * (im.0 + im.1);
        (re, (centre - target / 2.0, centre + target / 2.0))
    }
}

/// A generous first guess at the extent, only used to place the probe.
fn probe_scale(method: &Method) -> f64 {
    let scale = match &method.kind {
        MethodKind::RungeKutta(tableau) => {
            let function = StabilityFunction::from_tableau(tableau);
            let poles = function.poles();
            if !poles.is_empty() {
                // The poles are where the picture has its structure.
                3.0 * poles.iter().fold(1.0f64, |acc, p| acc.max(p.abs()))
            } else {
                // An explicit method is bounded by its real axis limit.
                let limit = function.real_stability_limit();
                if limit.is_finite() {
                    2.5 * limit.abs()
                } else {
                    12.0
                }
            }
        }
        MethodKind::LinearMultistep(family) => match family.uniform_coefficients() {
            Ok(coefficients) => {
                // The boundary locus traces the edge of the region directly, but
                // it runs to infinity whenever sigma has a root on the unit
                // circle, which is exactly the case for the methods whose region
                // is a half plane. A middle quantile is immune to that tail and
                // still tracks the size of the picture.
                let polynomials = GeneratingPolynomials::from_coefficients(&coefficients);
                let mut moduli: Vec<f64> = polynomials
                    .boundary_locus(240)
                    .iter()
                    .map(|z| z.abs())
                    .filter(|v| v.is_finite())
                    .collect();
                if moduli.is_empty() {
                    return 12.0;
                }
                moduli.sort_by(f64::total_cmp);
                let median = moduli[moduli.len() / 2];
                4.0 * median
            }
            Err(_) => 12.0,
        },
    };
    scale.clamp(2.0, 200.0)
}
