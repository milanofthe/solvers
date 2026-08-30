//! Nonlinear stability and monotonicity.
//!
//! Linear stability asks what a method does to `y' = lambda y`. It is the
//! question everything else in this crate answers, and it is not the only one.
//! Two properties of a tableau say what happens on a problem that is not
//! linear, both of them derivable and neither of them usually tabulated.
//!
//! **Algebraic stability.** For a problem with a one sided Lipschitz constant
//! at most zero, two solutions never move apart. A method preserves that when
//! `b >= 0` and
//!
//! ```text
//! M = diag(b) A + A^T diag(b) - b b^T
//! ```
//!
//! is positive semidefinite. Algebraic stability implies B-stability, so this
//! one eigenvalue problem on an `s` by `s` matrix settles the whole question of
//! nonlinear contractivity. Gauss, Radau IA and IIA and Lobatto IIIC are
//! algebraically stable; the diagonally implicit families mostly are not.
//!
//! **The radius of absolute monotonicity.** The SSP coefficient of a method is
//! the largest `r` for which the step can be written as a convex combination of
//! forward Euler steps of size `h/r`. That is exactly
//!
//! ```text
//! (I + rA)^-1 exists,   rA(I + rA)^-1 >= 0,   r b^T (I + rA)^-1 >= 0,
//! (I + rA)^-1 e >= 0
//! ```
//!
//! entrywise, and the conditions are monotone in `r`, so a bisection finds it.
//! A method with `C > 0` keeps whatever bound forward Euler keeps, provided the
//! step is at most `C` times the one forward Euler needs, which is what makes
//! the SSP families worth their extra stages.
//!
//! References
//! ----------
//! * K. Burrage, J. C. Butcher, "Stability criteria for implicit Runge-Kutta
//!   methods", SIAM J. Numer. Anal. 16, 1979, doi:10.1137/0716037
//! * S. Gottlieb, D. Ketcheson, C.-W. Shu, "Strong Stability Preserving
//!   Runge-Kutta and Multistep Time Discretizations", World Scientific 2011,
//!   doi:10.1142/7498

use crate::linalg::{symmetric_eigenvalues, Lu, Matrix};
use crate::method::RkTableau;

/// Whether the tableau is algebraically stable, and by how much.
///
/// The margin is the smallest eigenvalue of `M`. It is reported because zero is
/// the interesting value: the Gauss methods sit exactly on it, with `M = 0`,
/// and a method a hair below it is a method whose coefficients are rounded
/// rather than one that fails.
#[derive(Copy, Clone, Debug, serde::Serialize)]
pub struct AlgebraicStability {
    pub stable: bool,
    pub smallest_eigenvalue: f64,
    pub weights_non_negative: bool,
}

pub fn algebraic_stability(tableau: &RkTableau) -> AlgebraicStability {
    let s = tableau.stages;
    let a = tableau.a.map(|v| v.value());
    let b: Vec<f64> = tableau.b.iter().map(|v| v.value()).collect();
    let scale = b.iter().fold(0.0f64, |acc, v| acc.max(v.abs())).max(1.0);

    let mut m = Matrix::<f64>::zeros(s, s);
    for i in 0..s {
        for j in 0..s {
            m[(i, j)] = b[i] * a[(i, j)] + b[j] * a[(j, i)] - b[i] * b[j];
        }
    }

    let smallest = symmetric_eigenvalues(&m)
        .into_iter()
        .fold(f64::INFINITY, f64::min);
    let non_negative = b.iter().all(|v| *v >= -1e-12 * scale);
    AlgebraicStability {
        stable: non_negative && smallest >= -1e-10 * scale * scale,
        smallest_eigenvalue: smallest,
        weights_non_negative: non_negative,
    }
}

/// The strong stability preserving coefficient, by bisection on the radius of
/// absolute monotonicity.
///
/// Zero for a method that has none, which is most of them. The bisection is
/// sound because the conditions are monotone: a method absolutely monotone at
/// `r` is absolutely monotone at everything below it.
pub fn ssp_coefficient(tableau: &RkTableau) -> f64 {
    if !absolutely_monotone(tableau, 1e-9) {
        return 0.0;
    }
    let ceiling = 1024.0;
    let mut low = 1e-9;
    let mut high = 1e-9;
    // Find a radius it fails at before bisecting between the two. A method
    // still monotone at the ceiling is monotone at every radius: backward
    // Euler and the implicit midpoint rule are, and reporting a number there
    // would be reporting where the search stopped.
    for _ in 0..40 {
        let candidate = high * 2.0;
        if candidate > ceiling {
            return f64::INFINITY;
        }
        if !absolutely_monotone(tableau, candidate) {
            high = candidate;
            break;
        }
        low = candidate;
        high = candidate;
    }
    if low == high {
        return low;
    }
    for _ in 0..60 {
        let middle = 0.5 * (low + high);
        if absolutely_monotone(tableau, middle) {
            low = middle;
        } else {
            high = middle;
        }
    }
    low
}

/// The four sign conditions at one radius.
fn absolutely_monotone(tableau: &RkTableau, r: f64) -> bool {
    let s = tableau.stages;
    let a = tableau.a.map(|v| v.value());
    let b: Vec<f64> = tableau.b.iter().map(|v| v.value()).collect();

    let mut shifted = Matrix::<f64>::zeros(s, s);
    for i in 0..s {
        for j in 0..s {
            shifted[(i, j)] = r * a[(i, j)];
        }
        shifted[(i, i)] += 1.0;
    }
    let lu = Lu::factor(shifted);
    if lu.is_singular() {
        return false;
    }

    // Columns of `(I + rA)^-1`, which is all three of the remaining conditions.
    let mut inverse = Matrix::<f64>::zeros(s, s);
    for column in 0..s {
        let mut unit = vec![0.0; s];
        unit[column] = 1.0;
        if !lu.solve_in_place(&mut unit) {
            return false;
        }
        for row in 0..s {
            inverse[(row, column)] = unit[row];
        }
    }

    let tolerance = -1e-12;
    // `(I + rA)^-1 e >= 0`.
    for row in 0..s {
        let sum: f64 = (0..s).map(|j| inverse[(row, j)]).sum();
        if sum < tolerance {
            return false;
        }
    }
    // `r A (I + rA)^-1 >= 0`.
    for i in 0..s {
        for j in 0..s {
            let entry: f64 = (0..s).map(|k| r * a[(i, k)] * inverse[(k, j)]).sum();
            if entry < tolerance {
                return false;
            }
        }
    }
    // `r b^T (I + rA)^-1 >= 0`.
    for j in 0..s {
        let entry: f64 = (0..s).map(|i| r * b[i] * inverse[(i, j)]).sum();
        if entry < tolerance {
            return false;
        }
    }
    true
}
