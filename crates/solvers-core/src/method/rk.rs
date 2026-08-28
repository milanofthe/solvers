//! Runge-Kutta tableaux.
//!
//! One representation covers explicit, diagonally implicit and fully implicit
//! methods. The structure is derived from the coefficients rather than declared,
//! so a method file only has to state `A`, `b` and optionally `c`, the embedded
//! weights and a dense output.

use super::coeff_serde::CoeffValue;
use crate::linalg::Matrix;
use crate::num::{Coeff, Field};
use serde::{Deserialize, Serialize};

/// Raw tableau as it appears in a method file.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct RkTableauFile {
    /// Stage coupling matrix. Rows may be truncated at the last nonzero entry.
    pub a: Vec<Vec<CoeffValue>>,
    /// Weights of the propagating solution.
    pub b: Vec<CoeffValue>,
    /// Abscissae. Defaults to the row sums of `a`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub c: Option<Vec<CoeffValue>>,
    /// Weights of the embedded solution used for error estimation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub b_embedded: Option<Vec<CoeffValue>>,
    /// Dense output as `b_i(theta) = sum_j dense[i][j] * theta^(j+1)`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dense_output: Option<Vec<Vec<CoeffValue>>>,
}

/// How the stages couple, which decides how a step is computed.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Structure {
    /// Strictly lower triangular: stages evaluate in sequence.
    Explicit,
    /// Lower triangular: one implicit solve of size `n` per stage.
    DiagonallyImplicit,
    /// Dense: one coupled implicit solve of size `s * n`.
    FullyImplicit,
}

/// A validated Runge-Kutta tableau together with its derived properties.
#[derive(Clone, Debug)]
pub struct RkTableau {
    pub stages: usize,
    pub a: Matrix<Coeff>,
    pub b: Vec<Coeff>,
    pub c: Vec<Coeff>,
    pub b_embedded: Option<Vec<Coeff>>,
    pub dense_output: Option<Matrix<Coeff>>,
    pub structure: Structure,
    /// True when every nonzero diagonal entry is the same (SDIRK, ESDIRK).
    pub singly_diagonal: bool,
    /// True when the first stage is explicit (ESDIRK).
    pub explicit_first_stage: bool,
    /// True when the last row of `a` equals `b`, so the step ends on the
    /// solution point. Gives L-stability a chance and enables FSAL.
    pub stiffly_accurate: bool,
    /// The shared diagonal value for singly diagonal methods.
    pub gamma: Option<Coeff>,
}

#[derive(Debug)]
pub struct TableauError(pub String);

impl std::fmt::Display for TableauError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for TableauError {}

fn tol_zero(c: Coeff) -> bool {
    match c {
        Coeff::Exact(r) => r.is_zero(),
        Coeff::Approx(v) => v.abs() < 1e-14,
    }
}

fn tol_eq(a: Coeff, b: Coeff) -> bool {
    match (a, b) {
        (Coeff::Exact(x), Coeff::Exact(y)) => x == y,
        _ => (a.value() - b.value()).abs() < 1e-13 * a.value().abs().max(b.value().abs()).max(1.0),
    }
}

impl RkTableau {
    pub fn from_file(file: &RkTableauFile) -> Result<RkTableau, TableauError> {
        let s = file.b.len();
        if s == 0 {
            return Err(TableauError("tableau has no stages".into()));
        }
        if file.a.len() != s {
            return Err(TableauError(format!(
                "matrix a has {} rows but b has {} entries",
                file.a.len(),
                s
            )));
        }

        let mut a = Matrix::<Coeff>::zeros(s, s);
        for (i, row) in file.a.iter().enumerate() {
            if row.len() > s {
                return Err(TableauError(format!("row {i} of a has more than {s} entries")));
            }
            for (j, v) in row.iter().enumerate() {
                a[(i, j)] = v.get();
            }
        }

        let b: Vec<Coeff> = file.b.iter().map(|v| v.get()).collect();

        let c: Vec<Coeff> = match &file.c {
            Some(c) => {
                if c.len() != s {
                    return Err(TableauError(format!("c has {} entries, expected {s}", c.len())));
                }
                c.iter().map(|v| v.get()).collect()
            }
            None => (0..s)
                .map(|i| {
                    let mut acc = Coeff::zero();
                    for j in 0..s {
                        acc = acc + a[(i, j)];
                    }
                    acc
                })
                .collect(),
        };

        let b_embedded = match &file.b_embedded {
            Some(be) => {
                if be.len() != s {
                    return Err(TableauError(format!(
                        "b_embedded has {} entries, expected {s}",
                        be.len()
                    )));
                }
                Some(be.iter().map(|v| v.get()).collect::<Vec<_>>())
            }
            None => None,
        };

        let dense_output = match &file.dense_output {
            Some(rows) => {
                if rows.len() != s {
                    return Err(TableauError(format!(
                        "dense_output has {} rows, expected {s}",
                        rows.len()
                    )));
                }
                let deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
                let mut m = Matrix::<Coeff>::zeros(s, deg);
                for (i, row) in rows.iter().enumerate() {
                    for (j, v) in row.iter().enumerate() {
                        m[(i, j)] = v.get();
                    }
                }
                Some(m)
            }
            None => None,
        };

        // Derive the structure from the sparsity pattern.
        let mut strictly_lower = true;
        let mut lower = true;
        for i in 0..s {
            for j in 0..s {
                if !tol_zero(a[(i, j)]) {
                    if j > i {
                        lower = false;
                        strictly_lower = false;
                    } else if j == i {
                        strictly_lower = false;
                    }
                }
            }
        }
        let structure = if strictly_lower {
            Structure::Explicit
        } else if lower {
            Structure::DiagonallyImplicit
        } else {
            Structure::FullyImplicit
        };

        let mut gamma = None;
        let mut singly_diagonal = structure == Structure::DiagonallyImplicit;
        if singly_diagonal {
            for i in 0..s {
                let d = a[(i, i)];
                if tol_zero(d) {
                    continue;
                }
                match gamma {
                    None => gamma = Some(d),
                    Some(g) => {
                        if !tol_eq(g, d) {
                            singly_diagonal = false;
                            break;
                        }
                    }
                }
            }
        }

        let explicit_first_stage = (0..s).all(|j| tol_zero(a[(0, j)])) && structure != Structure::Explicit;
        let stiffly_accurate = (0..s).all(|j| tol_eq(a[(s - 1, j)], b[j]));

        Ok(RkTableau {
            stages: s,
            a,
            b,
            c,
            b_embedded,
            dense_output,
            structure,
            singly_diagonal,
            explicit_first_stage,
            stiffly_accurate,
            gamma: if singly_diagonal { gamma } else { None },
        })
    }

    pub fn is_explicit(&self) -> bool {
        self.structure == Structure::Explicit
    }

    pub fn has_embedded(&self) -> bool {
        self.b_embedded.is_some()
    }

    /// Explicit methods can reuse the last stage as the first of the next step
    /// when the last row of `a` equals `b` and the last abscissa is one.
    pub fn is_fsal(&self) -> bool {
        self.is_explicit()
            && self.stages > 1
            && self.stiffly_accurate
            && tol_eq(self.c[self.stages - 1], Coeff::one())
    }

    /// Error weights `b - b_embedded`.
    pub fn error_weights(&self) -> Option<Vec<Coeff>> {
        self.b_embedded
            .as_ref()
            .map(|be| self.b.iter().zip(be).map(|(x, y)| *x - *y).collect())
    }

    /// Float view used by the hot loop.
    pub fn runtime(&self) -> RkRuntime {
        RkRuntime {
            stages: self.stages,
            a: self.a.map(|v| v.value()),
            b: self.b.iter().map(|v| v.value()).collect(),
            c: self.c.iter().map(|v| v.value()).collect(),
            e: self
                .error_weights()
                .map(|e| e.iter().map(|v| v.value()).collect()),
            dense: self.dense_output.as_ref().map(|d| d.map(|v| v.value())),
            structure: self.structure,
            gamma: self.gamma.map(|g| g.value()),
            explicit_first_stage: self.explicit_first_stage,
            stiffly_accurate: self.stiffly_accurate,
            fsal: self.is_fsal(),
        }
    }
}

/// Tableau reduced to floats, ready for time stepping.
#[derive(Clone, Debug)]
pub struct RkRuntime {
    pub stages: usize,
    pub a: Matrix<f64>,
    pub b: Vec<f64>,
    pub c: Vec<f64>,
    /// Error weights `b - b_embedded`, if an embedded method exists.
    pub e: Option<Vec<f64>>,
    pub dense: Option<Matrix<f64>>,
    pub structure: Structure,
    pub gamma: Option<f64>,
    pub explicit_first_stage: bool,
    pub stiffly_accurate: bool,
    pub fsal: bool,
}
