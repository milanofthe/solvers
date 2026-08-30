//! Rosenbrock-Wanner tableaux.
//!
//! A ROW method is linearly implicit: the Jacobian is part of the formula
//! rather than a device for solving one, so a step costs one factorization and
//! `s` back substitutions and never iterates. Wanner's form is
//!
//! ```text
//! (I - h gamma_ii J) K_i = h f(t + c_i h, y + sum_j alpha_ij K_j)
//!                          + h J sum_j gamma_ij K_j
//!                          + h^2 d_i df/dt
//! y1 = y + sum_i b_i K_i
//! ```
//!
//! with `alpha` strictly lower triangular, `gamma` lower triangular, `J` the
//! Jacobian at the step's own starting point, and
//! `c_i = sum_j alpha_ij`, `d_i = sum_j gamma_ij`. Those two row sums are not
//! free: they are what makes the method behave on a non-autonomous problem the
//! way the autonomous order conditions say it should.
//!
//! Implementations do not run that form. They run the equivalent one in which
//! the matrix-vector product with `J` has been substituted away, which is what
//! the published coefficient tables are usually written in:
//!
//! ```text
//! (I/(h gamma) - J) S_i = f(t + c_i h, y + sum_j a_ij S_j)
//!                         + sum_j (C_ij / h) S_j + h d_i df/dt
//! y1 = y + sum_i m_i S_i
//! ```
//!
//! with `S = Gamma K` and therefore
//! `a = alpha Gamma^-1`, `C = diag(1/gamma_ii) - Gamma^-1`, `m^T = b^T Gamma^-1`.
//! The file stores `alpha`, `gamma` and `b`, which is the form the order
//! conditions are stated in; the running coefficients are derived here.
//!
//! Reference: J. Rang, "Rosenbrock-Wanner methods: construction and mission",
//! and E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II",
//! 2nd ed., Springer 1996, doi:10.1007/978-3-642-05221-7

use super::coeff_serde::CoeffValue;
use crate::linalg::Matrix;
use crate::num::{Coeff, Field};
use serde::{Deserialize, Serialize};

/// Raw Rosenbrock coefficients as they appear in a method file.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct RosenbrockFile {
    /// Stage coupling of the function argument, strictly lower triangular.
    /// Rows may be truncated at the last nonzero entry.
    pub alpha: Vec<Vec<CoeffValue>>,
    /// Coupling through the Jacobian, lower triangular including the diagonal.
    pub gamma: Vec<Vec<CoeffValue>>,
    /// Weights of the propagating solution.
    pub b: Vec<CoeffValue>,
    /// Weights of the embedded solution used for error estimation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub b_embedded: Option<Vec<CoeffValue>>,
}

/// A validated Rosenbrock tableau together with its derived properties.
#[derive(Clone, Debug)]
pub struct RosenbrockTableau {
    pub stages: usize,
    pub alpha: Matrix<Coeff>,
    pub gamma: Matrix<Coeff>,
    pub b: Vec<Coeff>,
    pub b_embedded: Option<Vec<Coeff>>,
    /// Abscissae `c_i = sum_j alpha_ij`.
    pub c: Vec<Coeff>,
    /// Time derivative weights `d_i = sum_j gamma_ij`.
    pub d: Vec<Coeff>,
    /// True when every diagonal entry of `gamma` is the same, which is what
    /// lets the whole step share one factorization.
    pub singly_diagonal: bool,
    /// The shared diagonal value, when there is one.
    pub diagonal: Option<Coeff>,
    /// True when `alpha_sj + gamma_sj = b_j` and the last abscissa is one, so
    /// the step ends on the solution point. This is what buys L-stability.
    pub stiffly_accurate: bool,
}

#[derive(Debug)]
pub struct RosenbrockError(pub String);

impl std::fmt::Display for RosenbrockError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for RosenbrockError {}

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

/// Read a triangular block into a square matrix, rejecting anything above the
/// diagonal the caller does not allow.
fn triangular(
    rows: &[Vec<CoeffValue>],
    s: usize,
    name: &str,
    strict: bool,
) -> Result<Matrix<Coeff>, RosenbrockError> {
    if rows.len() != s {
        return Err(RosenbrockError(format!(
            "{name} has {} rows but b has {s} entries",
            rows.len()
        )));
    }
    let mut m = Matrix::<Coeff>::zeros(s, s);
    for (i, row) in rows.iter().enumerate() {
        if row.len() > s {
            return Err(RosenbrockError(format!(
                "row {i} of {name} has more than {s} entries"
            )));
        }
        for (j, v) in row.iter().enumerate() {
            let value = v.get();
            let above = if strict { j >= i } else { j > i };
            if above && !tol_zero(value) {
                return Err(RosenbrockError(format!(
                    "{name}[{i}][{j}] is nonzero, which the structure does not allow"
                )));
            }
            m[(i, j)] = value;
        }
    }
    Ok(m)
}

impl RosenbrockTableau {
    pub fn from_file(file: &RosenbrockFile) -> Result<RosenbrockTableau, RosenbrockError> {
        let s = file.b.len();
        if s == 0 {
            return Err(RosenbrockError("tableau has no stages".into()));
        }

        let alpha = triangular(&file.alpha, s, "alpha", true)?;
        let gamma = triangular(&file.gamma, s, "gamma", false)?;
        for i in 0..s {
            if tol_zero(gamma[(i, i)]) {
                return Err(RosenbrockError(format!(
                    "gamma[{i}][{i}] is zero, so stage {i} has no linear system to solve"
                )));
            }
        }

        let b: Vec<Coeff> = file.b.iter().map(|v| v.get()).collect();
        let b_embedded = match &file.b_embedded {
            Some(be) if be.len() != s => {
                return Err(RosenbrockError(format!(
                    "b_embedded has {} entries, expected {s}",
                    be.len()
                )))
            }
            Some(be) => Some(be.iter().map(|v| v.get()).collect::<Vec<_>>()),
            None => None,
        };

        let row_sums = |m: &Matrix<Coeff>| -> Vec<Coeff> {
            (0..s)
                .map(|i| (0..s).fold(Coeff::zero(), |acc, j| acc + m[(i, j)]))
                .collect()
        };
        let c = row_sums(&alpha);
        let d = row_sums(&gamma);

        let diagonal = {
            let first = gamma[(0, 0)];
            (0..s)
                .all(|i| tol_eq(gamma[(i, i)], first))
                .then_some(first)
        };

        let stiffly_accurate = tol_eq(c[s - 1], Coeff::one())
            && (0..s).all(|j| tol_eq(alpha[(s - 1, j)] + gamma[(s - 1, j)], b[j]));

        Ok(RosenbrockTableau {
            stages: s,
            alpha,
            gamma,
            b,
            b_embedded,
            c,
            d,
            singly_diagonal: diagonal.is_some(),
            diagonal,
            stiffly_accurate,
        })
    }

    pub fn has_embedded(&self) -> bool {
        self.b_embedded.is_some()
    }

    /// `alpha + gamma`, the coefficient matrix of the diagonally implicit
    /// Runge-Kutta method this one is a single Newton step of. It carries the
    /// whole linear stability of the method, which is why a ROW method has the
    /// stability function of a DIRK.
    pub fn implicit_matrix(&self) -> Matrix<Coeff> {
        let s = self.stages;
        let mut m = Matrix::<Coeff>::zeros(s, s);
        for i in 0..s {
            for j in 0..s {
                m[(i, j)] = self.alpha[(i, j)] + self.gamma[(i, j)];
            }
        }
        m
    }

    /// Float view in the substituted coordinates the stepper runs in.
    pub fn runtime(&self) -> RosenbrockRuntime {
        let s = self.stages;
        let gamma = self.gamma.map(|v| v.value());
        let inverse = invert_lower(&gamma);

        let alpha = self.alpha.map(|v| v.value());
        let mut a = Matrix::<f64>::zeros(s, s);
        for i in 0..s {
            for j in 0..s {
                a[(i, j)] = (0..s).map(|k| alpha[(i, k)] * inverse[(k, j)]).sum();
            }
        }

        let mut c_matrix = Matrix::<f64>::zeros(s, s);
        for i in 0..s {
            for j in 0..s {
                c_matrix[(i, j)] = -inverse[(i, j)];
            }
            c_matrix[(i, i)] += 1.0 / gamma[(i, i)];
        }

        let weights = |w: &[Coeff]| -> Vec<f64> {
            (0..s)
                .map(|j| (0..s).map(|k| w[k].value() * inverse[(k, j)]).sum())
                .collect()
        };
        let m = weights(&self.b);
        let e = self
            .b_embedded
            .as_ref()
            .map(|be| weights(be))
            .map(|me| m.iter().zip(&me).map(|(x, y)| x - y).collect());

        RosenbrockRuntime {
            stages: s,
            a,
            c_matrix,
            m,
            e,
            c: self.c.iter().map(|v| v.value()).collect(),
            d: self.d.iter().map(|v| v.value()).collect(),
            diagonal: self.diagonal.map(|v| v.value()),
            stiffly_accurate: self.stiffly_accurate,
        }
    }
}

/// Invert a lower triangular matrix by forward substitution, column by column.
fn invert_lower(m: &Matrix<f64>) -> Matrix<f64> {
    let n = m.rows();
    let mut inverse = Matrix::<f64>::zeros(n, n);
    for column in 0..n {
        for row in column..n {
            let mut acc = if row == column { 1.0 } else { 0.0 };
            for k in column..row {
                acc -= m[(row, k)] * inverse[(k, column)];
            }
            inverse[(row, column)] = acc / m[(row, row)];
        }
    }
    inverse
}

/// Coefficients reduced to floats, in the substituted form the stepper uses.
#[derive(Clone, Debug)]
pub struct RosenbrockRuntime {
    pub stages: usize,
    /// `alpha Gamma^-1`, the coupling of the function argument.
    pub a: Matrix<f64>,
    /// `diag(1/gamma_ii) - Gamma^-1`, the coupling of the past stages.
    pub c_matrix: Matrix<f64>,
    /// `b^T Gamma^-1`, the weights of the propagating solution.
    pub m: Vec<f64>,
    /// `m - m_embedded`, the error weights, when an embedded pair exists.
    pub e: Option<Vec<f64>>,
    /// Abscissae.
    pub c: Vec<f64>,
    /// Time derivative weights.
    pub d: Vec<f64>,
    pub diagonal: Option<f64>,
    pub stiffly_accurate: bool,
}

impl RosenbrockRuntime {
    /// Whether one factorization serves the whole step.
    pub fn shares_one_factorization(&self) -> bool {
        self.diagonal.is_some()
    }
}
