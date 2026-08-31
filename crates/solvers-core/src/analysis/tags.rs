//! Tags derived from what a method actually is.
//!
//! The library is browsed by property, not by name, so the properties have to
//! be first class and they have to be derived. A tag here is never a label
//! somebody typed into a file; it follows from the coefficients, which means a
//! filter on "L-stable" returns exactly the methods that are.

use super::MethodReport;
use crate::method::{Method, MethodKind};

/// A tag with the group it belongs to, so the interface can offer them in
/// coherent sets rather than as one flat list.
#[derive(Clone, Debug, serde::Serialize)]
pub struct Tag {
    pub group: &'static str,
    pub name: String,
}

fn tag(group: &'static str, name: &str) -> Tag {
    Tag {
        group,
        name: name.to_string(),
    }
}

/// Every tag that applies to a method.
pub fn tags(method: &Method, report: &MethodReport) -> Vec<Tag> {
    let mut tags = Vec::new();

    // Structure: how a step is actually computed.
    match &method.kind {
        MethodKind::Additive(pair) => {
            tags.push(tag("structure", "additive"));
            if pair.implicit.singly_diagonal {
                tags.push(tag("structure", "singly diagonal"));
            }
            if pair.implicit.explicit_first_stage {
                tags.push(tag("structure", "explicit first stage"));
            }
            if pair.implicit.stiffly_accurate {
                tags.push(tag("structure", "stiffly accurate"));
            }
        }
        MethodKind::RungeKutta(tableau) => {
            if tableau.is_explicit() {
                tags.push(tag("structure", "explicit"));
            } else if tableau.singly_diagonal {
                tags.push(tag("structure", "singly diagonal"));
            } else if tableau.structure == crate::method::Structure::DiagonallyImplicit {
                tags.push(tag("structure", "diagonally implicit"));
            } else {
                tags.push(tag("structure", "fully implicit"));
            }
            if tableau.explicit_first_stage {
                tags.push(tag("structure", "explicit first stage"));
            }
            if tableau.is_fsal() {
                tags.push(tag("structure", "FSAL"));
            }
            if tableau.stiffly_accurate {
                tags.push(tag("structure", "stiffly accurate"));
            }
            if tableau.has_embedded() {
                tags.push(tag("control", "embedded pair"));
            } else {
                tags.push(tag("control", "step doubling"));
            }
        }
        MethodKind::LinearMultistep(family) => {
            tags.push(tag("structure", "multistep"));
            if family.implicit {
                tags.push(tag("structure", "implicit"));
            } else {
                tags.push(tag("structure", "explicit"));
            }
            if family.startup.is_some() {
                tags.push(tag("structure", "needs a starter"));
            } else {
                tags.push(tag("structure", "self starting"));
            }
            tags.push(tag("control", "order reduced estimate"));
        }
        MethodKind::Rosenbrock(tableau) => {
            tags.push(tag("structure", "linearly implicit"));
            if tableau.singly_diagonal {
                tags.push(tag("structure", "singly diagonal"));
            }
            if tableau.stiffly_accurate {
                tags.push(tag("structure", "stiffly accurate"));
            }
            // The Jacobian is the formula here, not an accelerator, so it is
            // worth saying out loud: this is the one family that cannot be run
            // without one.
            tags.push(tag("structure", "needs a Jacobian"));
            if tableau.has_embedded() {
                tags.push(tag("control", "embedded pair"));
            } else {
                tags.push(tag("control", "step doubling"));
            }
        }
    }

    // Stability: the property that decides whether a method can be used at all
    // on a given problem.
    if report.l_stable {
        tags.push(tag("stability", "L-stable"));
    } else if report.a_stable {
        tags.push(tag("stability", "A-stable"));
    } else if let Some(angle) = report.alpha_angle {
        if angle >= 1.0 {
            tags.push(tag("stability", "A(alpha)-stable"));
        } else {
            tags.push(tag("stability", "conditionally stable"));
        }
    } else {
        tags.push(tag("stability", "conditionally stable"));
    }
    if report.a_stable || report.alpha_angle.map_or(false, |a| a >= 45.0) {
        tags.push(tag("use", "stiff problems"));
    } else {
        tags.push(tag("use", "non stiff problems"));
    }

    // Nonlinear stability, which linear stability says nothing about.
    if report.algebraically_stable == Some(true) {
        tags.push(tag("stability", "algebraically stable"));
    }
    let preserves_strong_stability = matches!(
        report.ssp_coefficient,
        Some(super::Limit::Unbounded(_)) | Some(super::Limit::Finite(1e-6..))
    );
    if preserves_strong_stability {
        tags.push(tag("stability", "strong stability preserving"));
    }
    if report.dissipation_order.is_none() && report.dispersion_order.is_some() {
        tags.push(tag("geometry", "non dissipative"));
    }

    // Geometric properties, which the analysis does not derive and the file
    // therefore has to claim.
    if method.properties.symplectic == Some(true) {
        tags.push(tag("geometry", "symplectic"));
    }
    if method.properties.symmetric == Some(true) {
        tags.push(tag("geometry", "symmetric"));
    }

    // Order, in the bands one actually chooses between.
    tags.push(tag(
        "order",
        match report.computed_order {
            0..=2 => "order 1 to 2",
            3..=4 => "order 3 to 4",
            5..=6 => "order 5 to 6",
            _ => "order 7 and above",
        },
    ));
    if report.adaptive {
        tags.push(tag("control", "adaptive"));
    } else {
        tags.push(tag("control", "fixed step"));
    }

    if report.exact_arithmetic {
        tags.push(tag("provenance", "exact coefficients"));
    }
    if method.references.iter().any(|r| r.doi.is_some()) {
        tags.push(tag("provenance", "published source"));
    }

    for extra in &method.tags {
        tags.push(tag("other", extra));
    }
    tags
}
