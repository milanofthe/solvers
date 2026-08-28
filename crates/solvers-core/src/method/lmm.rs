//! Linear multistep methods.
//!
//! A k-step method is written as
//!
//! ```text
//! sum_{j=0..k} alpha_j * y_{n-j} = h_n * sum_{j=0..k} beta_j * f_{n-j}
//! ```
//!
//! and a method file only declares which of the `alpha_j` and `beta_j` are
//! fixed and which are free. The free ones are then determined from the
//! exactness conditions on the actual step size history, which makes the
//! variable step form of BDF, Adams-Bashforth, Adams-Moulton, Nystrom and
//! Milne-Simpson fall out of a single implementation.

use super::coeff_serde::Slot;
use crate::linalg::{Lu, Matrix};
use serde::{Deserialize, Serialize};

/// Which coefficient is pinned to break the scaling degeneracy of a
/// homogeneous coefficient system.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Deserialize, Serialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Normalization {
    /// Rescale the solved coefficients so that `alpha_0 = 1`.
    #[default]
    Alpha0,
    /// Leave the coefficients as solved.
    None,
}

/// Raw multistep description from a method file.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct LmmFile {
    /// Number of steps `k`.
    pub steps: usize,
    /// Coefficients of `y_{n-j}`, index `j = 0..k`.
    pub alpha: Vec<Slot>,
    /// Coefficients of `f_{n-j}`, index `j = 0..k`.
    pub beta: Vec<Slot>,
    #[serde(default)]
    pub normalization: Normalization,
    /// Method id used to generate the missing history at start up.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub startup: Option<String>,
    /// Lowest step count the family is still a valid method at.
    ///
    /// One for BDF and Adams, which have a member at every step count down to
    /// one. Two for Nystrom and Milne-Simpson, whose alpha pattern stops making
    /// sense when truncated. It bounds both the variable order driver and the
    /// order reduced error estimate.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_steps: Option<usize>,
}

/// A validated multistep family, parametrized by its free coefficient pattern.
#[derive(Clone, Debug)]
pub struct LmmFamily {
    pub steps: usize,
    pub alpha: Vec<Slot>,
    pub beta: Vec<Slot>,
    pub normalization: Normalization,
    pub startup: Option<String>,
    pub min_steps: usize,
    /// Whether `beta_0` can be nonzero, i.e. whether a step needs a solve.
    pub implicit: bool,
}

#[derive(Debug)]
pub struct LmmError(pub String);

impl std::fmt::Display for LmmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for LmmError {}

/// Coefficients for one concrete step size history.
#[derive(Clone, Debug, PartialEq)]
pub struct LmmCoefficients {
    pub alpha: Vec<f64>,
    pub beta: Vec<f64>,
    /// Highest degree of polynomial the formula integrates exactly.
    pub order: usize,
}

impl LmmCoefficients {
    pub fn is_implicit(&self) -> bool {
        self.beta[0] != 0.0
    }
}

/// `theta^d` with the `0^0 = 1` convention the exactness conditions need.
fn pow(theta: f64, d: i32) -> f64 {
    if d < 0 {
        return 0.0;
    }
    if d == 0 {
        1.0
    } else {
        theta.powi(d)
    }
}

impl LmmFamily {
    pub fn from_file(file: &LmmFile) -> Result<LmmFamily, LmmError> {
        let k = file.steps;
        if k == 0 {
            return Err(LmmError("a multistep method needs at least one step".into()));
        }
        if file.alpha.len() != k + 1 || file.beta.len() != k + 1 {
            return Err(LmmError(format!(
                "alpha and beta must have {} entries for a {}-step method",
                k + 1,
                k
            )));
        }
        let implicit = !matches!(file.beta[0], Slot::Fixed(c) if c.value() == 0.0);
        Ok(LmmFamily {
            steps: k,
            alpha: file.alpha.clone(),
            beta: file.beta.clone(),
            normalization: file.normalization,
            startup: file.startup.clone(),
            min_steps: file.min_steps.unwrap_or(1).max(1),
            implicit,
        })
    }

    /// The same family with one step fewer, used for the error estimate and
    /// for the start up phase where the history is still short.
    pub fn with_steps(&self, k: usize) -> Option<LmmFamily> {
        if k == 0 || k > self.steps {
            return None;
        }
        Some(LmmFamily {
            steps: k,
            alpha: self.alpha[..=k].to_vec(),
            beta: self.beta[..=k].to_vec(),
            normalization: self.normalization,
            startup: self.startup.clone(),
            min_steps: self.min_steps,
            implicit: self.implicit,
        })
    }

    /// Normalized time offsets of the history points relative to the current
    /// step, `theta_j = (t_{n-j} - t_n) / h_n`.
    ///
    /// `steps[0]` is the current step size, `steps[i]` the one `i` steps back.
    pub fn thetas(&self, steps: &[f64]) -> Vec<f64> {
        let h = steps[0];
        let mut theta = vec![0.0; self.steps + 1];
        let mut acc = 0.0;
        for j in 1..=self.steps {
            acc += steps[(j - 1).min(steps.len() - 1)];
            theta[j] = -acc / h;
        }
        theta
    }

    /// Solve the exactness conditions for the free coefficients.
    ///
    /// Conditions are imposed for increasing polynomial degree; a degree whose
    /// row is empty in the unknowns is skipped but checked for consistency, so
    /// a malformed family shows up as an error rather than as a wrong order.
    pub fn coefficients(&self, steps: &[f64]) -> Result<LmmCoefficients, LmmError> {
        let k = self.steps;
        let theta = self.thetas(steps);

        // Index the unknowns: alphas first, then betas.
        let mut unknowns: Vec<(bool, usize)> = Vec::new();
        for (j, slot) in self.alpha.iter().enumerate() {
            if slot.is_free() {
                unknowns.push((true, j));
            }
        }
        for (j, slot) in self.beta.iter().enumerate() {
            if slot.is_free() {
                unknowns.push((false, j));
            }
        }
        let m = unknowns.len();
        if m == 0 {
            let alpha: Vec<f64> = self.alpha.iter().map(|s| s.fixed_value()).collect();
            let beta: Vec<f64> = self.beta.iter().map(|s| s.fixed_value()).collect();
            let order = self.attained_order(&alpha, &beta, &theta);
            return Ok(self.normalized(alpha, beta, order));
        }

        let mut rows: Vec<Vec<f64>> = Vec::with_capacity(m);
        let mut rhs: Vec<f64> = Vec::with_capacity(m);

        let mut degree = 0i32;
        let max_degree = (2 * (k + 2)) as i32;
        while rows.len() < m && degree <= max_degree {
            let mut row = vec![0.0; m];
            for (u, &(is_alpha, j)) in unknowns.iter().enumerate() {
                row[u] = if is_alpha {
                    pow(theta[j], degree)
                } else {
                    -(degree as f64) * pow(theta[j], degree - 1)
                };
            }
            // Move the fixed contributions to the right hand side.
            let mut fixed = 0.0;
            for (j, slot) in self.alpha.iter().enumerate() {
                if let Slot::Fixed(c) = slot {
                    fixed += c.value() * pow(theta[j], degree);
                }
            }
            for (j, slot) in self.beta.iter().enumerate() {
                if let Slot::Fixed(c) = slot {
                    fixed -= c.value() * (degree as f64) * pow(theta[j], degree - 1);
                }
            }

            let scale = row.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
            if scale < 1e-14 {
                // Degenerate row: no unknown appears. It must already hold.
                if fixed.abs() > 1e-9 {
                    return Err(LmmError(format!(
                        "inconsistent family: degree {degree} condition cannot be satisfied"
                    )));
                }
            } else {
                rows.push(row);
                rhs.push(-fixed);
            }
            degree += 1;
        }

        if rows.len() < m {
            return Err(LmmError("not enough independent order conditions".into()));
        }

        let matrix = Matrix::from_rows(&rows);
        let solution = Lu::factor(matrix)
            .solve(&rhs)
            .ok_or_else(|| LmmError("singular coefficient system".into()))?;

        let mut alpha: Vec<f64> = self.alpha.iter().map(|s| s.fixed_value()).collect();
        let mut beta: Vec<f64> = self.beta.iter().map(|s| s.fixed_value()).collect();
        for (u, &(is_alpha, j)) in unknowns.iter().enumerate() {
            if is_alpha {
                alpha[j] = solution[u];
            } else {
                beta[j] = solution[u];
            }
        }

        let order = self.attained_order(&alpha, &beta, &theta);
        Ok(self.normalized(alpha, beta, order))
    }

    fn normalized(&self, mut alpha: Vec<f64>, mut beta: Vec<f64>, order: usize) -> LmmCoefficients {
        if self.normalization == Normalization::Alpha0 && alpha[0].abs() > 1e-300 {
            let s = alpha[0];
            for v in alpha.iter_mut() {
                *v /= s;
            }
            for v in beta.iter_mut() {
                *v /= s;
            }
        }
        LmmCoefficients { alpha, beta, order }
    }

    /// Largest degree up to which the exactness conditions hold.
    fn attained_order(&self, alpha: &[f64], beta: &[f64], theta: &[f64]) -> usize {
        let magnitude = alpha
            .iter()
            .chain(beta.iter())
            .fold(1.0f64, |acc, v| acc.max(v.abs()));
        let mut order = 0usize;
        for degree in 0..=(2 * self.steps + 4) as i32 {
            let mut residual = 0.0;
            for (j, a) in alpha.iter().enumerate() {
                residual += a * pow(theta[j], degree);
            }
            for (j, b) in beta.iter().enumerate() {
                residual -= b * (degree as f64) * pow(theta[j], degree - 1);
            }
            if residual.abs() > 1e-9 * magnitude {
                return order;
            }
            order = degree.max(0) as usize;
        }
        order
    }

    /// Coefficients on a uniform grid, which is what the classical published
    /// tables list and what the stability analysis uses.
    pub fn uniform_coefficients(&self) -> Result<LmmCoefficients, LmmError> {
        let steps = vec![1.0; self.steps + 1];
        self.coefficients(&steps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::num::Coeff;

    fn family(k: usize, alpha: Vec<Slot>, beta: Vec<Slot>) -> LmmFamily {
        LmmFamily {
            steps: k,
            alpha,
            beta,
            normalization: Normalization::Alpha0,
            startup: None,
            min_steps: 1,
            implicit: true,
        }
    }

    fn fixed(v: f64) -> Slot {
        Slot::Fixed(Coeff::from_f64_rationalized(v))
    }

    #[test]
    fn bdf2_matches_the_textbook() {
        // All alphas free, beta_0 pinned to 1 before normalization.
        let f = family(
            2,
            vec![Slot::Free, Slot::Free, Slot::Free],
            vec![fixed(1.0), fixed(0.0), fixed(0.0)],
        );
        let c = f.uniform_coefficients().unwrap();
        // y_n - 4/3 y_{n-1} + 1/3 y_{n-2} = 2/3 h f_n
        assert!((c.alpha[0] - 1.0).abs() < 1e-12);
        assert!((c.alpha[1] + 4.0 / 3.0).abs() < 1e-12);
        assert!((c.alpha[2] - 1.0 / 3.0).abs() < 1e-12);
        assert!((c.beta[0] - 2.0 / 3.0).abs() < 1e-12);
        assert_eq!(c.order, 2);
    }

    #[test]
    fn adams_bashforth_3_matches_the_textbook() {
        let f = family(
            3,
            vec![fixed(1.0), fixed(-1.0), fixed(0.0), fixed(0.0)],
            vec![fixed(0.0), Slot::Free, Slot::Free, Slot::Free],
        );
        let c = f.uniform_coefficients().unwrap();
        assert!((c.beta[1] - 23.0 / 12.0).abs() < 1e-12);
        assert!((c.beta[2] + 16.0 / 12.0).abs() < 1e-12);
        assert!((c.beta[3] - 5.0 / 12.0).abs() < 1e-12);
        assert_eq!(c.order, 3);
    }

    #[test]
    fn adams_moulton_3_matches_the_textbook() {
        let f = family(
            3,
            vec![fixed(1.0), fixed(-1.0), fixed(0.0), fixed(0.0)],
            vec![Slot::Free, Slot::Free, Slot::Free, Slot::Free],
        );
        let c = f.uniform_coefficients().unwrap();
        // 9/24, 19/24, -5/24, 1/24
        assert!((c.beta[0] - 9.0 / 24.0).abs() < 1e-12);
        assert!((c.beta[1] - 19.0 / 24.0).abs() < 1e-12);
        assert!((c.beta[2] + 5.0 / 24.0).abs() < 1e-12);
        assert!((c.beta[3] - 1.0 / 24.0).abs() < 1e-12);
        assert_eq!(c.order, 4);
    }

    #[test]
    fn variable_step_bdf2_reduces_to_uniform() {
        let f = family(
            2,
            vec![Slot::Free, Slot::Free, Slot::Free],
            vec![fixed(1.0), fixed(0.0), fixed(0.0)],
        );
        let uniform = f.uniform_coefficients().unwrap();
        let same = f.coefficients(&[0.25, 0.25, 0.25]).unwrap();
        for (a, b) in uniform.alpha.iter().zip(&same.alpha) {
            assert!((a - b).abs() < 1e-12);
        }
        // A stretched history must still be second order exact.
        let stretched = f.coefficients(&[0.1, 0.3, 0.7]).unwrap();
        assert_eq!(stretched.order, 2);
    }

    #[test]
    fn backward_euler_is_the_one_step_case() {
        let f = family(1, vec![Slot::Free, Slot::Free], vec![fixed(1.0), fixed(0.0)]);
        let c = f.uniform_coefficients().unwrap();
        assert!((c.alpha[0] - 1.0).abs() < 1e-12);
        assert!((c.alpha[1] + 1.0).abs() < 1e-12);
        assert!((c.beta[0] - 1.0).abs() < 1e-12);
        assert_eq!(c.order, 1);
    }
}
