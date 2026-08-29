//! Linear stability.
//!
//! For a Runge-Kutta method the stability function is the rational function
//!
//! ```text
//! R(z) = det(I - zA + z e b^T) / det(I - zA)
//! ```
//!
//! Both determinants are characteristic polynomials in disguise, so they are
//! obtained exactly from the tableau with Faddeev-LeVerrier rather than by
//! evaluating a determinant on a grid. That gives the poles, the value at
//! infinity and an exact A-stability test through the E-polynomial, instead of
//! a picture that has to be eyeballed.
//!
//! For a linear multistep method the same questions are answered through the
//! root condition on `rho(w) - z sigma(w)`.
//!
//! References
//! ----------
//! * E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II",
//!   2nd ed., Springer 1996, doi:10.1007/978-3-642-05221-7
//! * G. Dahlquist, "A special stability problem for linear multistep methods",
//!   BIT 3, 1963, doi:10.1007/BF01963532

use crate::linalg::{poly_roots, Matrix};
use crate::method::{LmmCoefficients, RkTableau};
use crate::num::Complex;
use serde::Serialize;

/// Coefficients of `det(lambda I - A)`, index `i` multiplies `lambda^i`.
///
/// Faddeev-LeVerrier, which needs only matrix products and traces and is exact
/// enough at the sizes a tableau has.
pub fn characteristic_polynomial(a: &Matrix<f64>) -> Vec<f64> {
    let n = a.rows();
    let mut c = vec![0.0; n + 1];
    c[n] = 1.0;
    if n == 0 {
        return c;
    }

    let mut m = Matrix::<f64>::identity(n);
    for k in 1..=n {
        if k > 1 {
            // M_k = A * M_{k-1} + c_{n-k+1} * I
            let mut next = Matrix::<f64>::zeros(n, n);
            for i in 0..n {
                for j in 0..n {
                    let mut acc = 0.0;
                    for p in 0..n {
                        acc += a[(i, p)] * m[(p, j)];
                    }
                    next[(i, j)] = acc;
                }
                next[(i, i)] += c[n - k + 1];
            }
            m = next;
        }
        let mut trace = 0.0;
        for i in 0..n {
            let mut acc = 0.0;
            for p in 0..n {
                acc += a[(i, p)] * m[(p, i)];
            }
            trace += acc;
        }
        c[n - k] = -trace / k as f64;
    }
    c
}

/// A rational stability function `numerator / denominator`, both in powers of
/// `z` with index `i` multiplying `z^i`.
#[derive(Clone, Debug, Serialize)]
pub struct StabilityFunction {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
}

impl StabilityFunction {
    /// Build `R(z)` from a Butcher tableau.
    pub fn from_tableau(tableau: &RkTableau) -> StabilityFunction {
        let s = tableau.stages;
        let a = tableau.a.map(|v| v.value());
        let b: Vec<f64> = tableau.b.iter().map(|v| v.value()).collect();

        // det(I - zA) = z^s * charpoly_A(1/z), so the coefficients just reverse.
        let ca = characteristic_polynomial(&a);
        let denominator: Vec<f64> = (0..=s).map(|j| ca[s - j]).collect();

        // The numerator is the same determinant for A - e b^T.
        let mut shifted = a.clone();
        for i in 0..s {
            for j in 0..s {
                shifted[(i, j)] -= b[j];
            }
        }
        let cb = characteristic_polynomial(&shifted);
        let numerator: Vec<f64> = (0..=s).map(|j| cb[s - j]).collect();

        StabilityFunction {
            numerator,
            denominator,
        }
    }

    pub fn eval(&self, z: Complex) -> Complex {
        let n = eval_poly(&self.numerator, z);
        let d = eval_poly(&self.denominator, z);
        n / d
    }

    /// `R(infinity)`, the damping of infinitely stiff modes.
    pub fn at_infinity(&self) -> f64 {
        let dn = degree(&self.numerator);
        let dd = degree(&self.denominator);
        if dn < dd {
            0.0
        } else if dn > dd {
            f64::INFINITY
        } else {
            self.numerator[dn] / self.denominator[dd]
        }
    }

    /// Poles of `R`, i.e. the roots of the denominator.
    pub fn poles(&self) -> Vec<Complex> {
        if degree(&self.denominator) == 0 {
            return Vec::new();
        }
        let coefficients: Vec<Complex> = self
            .denominator
            .iter()
            .map(|v| Complex::real(*v))
            .collect();
        poly_roots(&coefficients)
    }

    /// Taylor coefficients of `R` around zero, up to `count` terms.
    pub fn taylor(&self, count: usize) -> Vec<f64> {
        // Power series division of the numerator by the denominator.
        let mut out = vec![0.0; count + 1];
        for k in 0..=count {
            let mut acc = *self.numerator.get(k).unwrap_or(&0.0);
            for j in 1..=k {
                acc -= self.denominator.get(j).copied().unwrap_or(0.0) * out[k - j];
            }
            out[k] = acc / self.denominator[0];
        }
        out
    }

    /// Largest `p` with `R(z) = exp(z) + O(z^(p+1))`.
    pub fn order_of_consistency(&self, max_order: usize) -> usize {
        let taylor = self.taylor(max_order + 1);
        let mut factorial = 1.0;
        let mut order = 0;
        for k in 0..=max_order {
            if k > 0 {
                factorial *= k as f64;
            }
            let target = 1.0 / factorial;
            if (taylor[k] - target).abs() > 1e-10 * target.abs().max(1.0) {
                return order;
            }
            order = k;
        }
        order
    }

    /// `|D(iy)|^2 - |N(iy)|^2` as a polynomial in `u = y^2`.
    ///
    /// A-stability on the imaginary axis is exactly the statement that this
    /// polynomial is non negative for all `u >= 0`.
    pub fn e_polynomial(&self) -> Vec<f64> {
        let d = squared_modulus_on_imaginary_axis(&self.denominator);
        let n = squared_modulus_on_imaginary_axis(&self.numerator);
        let len = d.len().max(n.len());
        (0..len)
            .map(|i| d.get(i).copied().unwrap_or(0.0) - n.get(i).copied().unwrap_or(0.0))
            .collect()
    }

    /// A-stable: no pole in the closed left half plane and `|R| <= 1` on the
    /// imaginary axis. The maximum modulus principle then gives `|R| <= 1`
    /// throughout the left half plane.
    pub fn is_a_stable(&self) -> bool {
        for pole in self.poles() {
            if pole.re <= 1e-10 {
                return false;
            }
        }
        let infinity = self.at_infinity();
        if !(infinity.abs() <= 1.0 + 1e-10) {
            return false;
        }
        polynomial_is_non_negative(&self.e_polynomial())
    }

    /// L-stable: A-stable and infinitely stiff modes are killed outright.
    pub fn is_l_stable(&self) -> bool {
        self.is_a_stable() && self.at_infinity().abs() <= 1e-12
    }

    /// Leftmost point of the real stability interval `[x, 0]`.
    pub fn real_stability_limit(&self) -> f64 {
        scan_limit(|x| self.eval(Complex::real(x)).abs() <= 1.0 + 1e-12, -1.0)
    }

    /// Largest `y` with the whole segment `[0, y]` of the imaginary axis stable.
    pub fn imaginary_stability_limit(&self) -> f64 {
        scan_limit(|y| self.eval(Complex::new(0.0, y)).abs() <= 1.0 + 1e-12, 1.0)
    }
}

fn degree(p: &[f64]) -> usize {
    let mut d = 0;
    for (i, v) in p.iter().enumerate() {
        if v.abs() > 1e-14 {
            d = i;
        }
    }
    d
}

fn eval_poly(p: &[f64], z: Complex) -> Complex {
    let mut acc = Complex::new(0.0, 0.0);
    for coefficient in p.iter().rev() {
        acc = acc * z + Complex::real(*coefficient);
    }
    acc
}

/// `|P(iy)|^2` as a polynomial in `u = y^2`, for a real coefficient `P`.
fn squared_modulus_on_imaginary_axis(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut in_y = vec![0.0; 2 * n];
    for j in 0..n {
        for k in 0..n {
            // i^(j-k) is real whenever j + k is even, which is the only case
            // that contributes.
            if (j + k) % 2 != 0 {
                continue;
            }
            let exponent = (j as i32 - k as i32).rem_euclid(4);
            let sign = if exponent == 0 { 1.0 } else { -1.0 };
            in_y[j + k] += sign * p[j] * p[k];
        }
    }
    (0..n).map(|m| in_y[2 * m]).collect()
}

/// Whether a real polynomial is non negative on `[0, infinity)`.
///
/// The sign can only change at a real root, so evaluating between consecutive
/// positive roots and beyond the largest one decides it.
fn polynomial_is_non_negative(p: &[f64]) -> bool {
    let trimmed: Vec<f64> = {
        let mut t = p.to_vec();
        while t.len() > 1 && t.last().map_or(false, |v| v.abs() < 1e-13) {
            t.pop();
        }
        t
    };
    if trimmed.len() == 1 {
        return trimmed[0] >= -1e-12;
    }

    let coefficients: Vec<Complex> = trimmed.iter().map(|v| Complex::real(*v)).collect();
    let mut breakpoints: Vec<f64> = poly_roots(&coefficients)
        .into_iter()
        .filter(|r| r.im.abs() < 1e-7 && r.re > 0.0)
        .map(|r| r.re)
        .collect();
    breakpoints.sort_by(|a, b| a.total_cmp(b));

    let mut samples = vec![0.0];
    let mut previous = 0.0;
    for point in &breakpoints {
        samples.push(0.5 * (previous + point));
        previous = *point;
    }
    samples.push(previous + 1.0);
    samples.push(previous * 2.0 + 10.0);

    let magnitude = trimmed.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
    samples.into_iter().all(|u| {
        let mut acc = 0.0;
        for coefficient in trimmed.iter().rev() {
            acc = acc * u + coefficient;
        }
        acc >= -1e-9 * magnitude
    })
}

/// Walk outwards from the origin along a direction until stability is lost.
///
/// Returns an infinite limit when the whole ray stays stable, which is the
/// honest answer for an A-stable method rather than wherever the scan stopped.
fn scan_limit(stable: impl Fn(f64) -> bool, direction: f64) -> f64 {
    let mut low = 0.0f64;
    let mut high = direction * 1e-6;
    // Grow until the first unstable point is found.
    for _ in 0..200 {
        if !stable(high) {
            break;
        }
        low = high;
        high *= 1.6;
        if high.abs() > 1e7 {
            return direction * f64::INFINITY;
        }
    }
    if stable(high) {
        return direction * f64::INFINITY;
    }
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if stable(mid) {
            low = mid;
        } else {
            high = mid;
        }
    }
    low
}

/// A sampled stability region, ready for plotting.
#[derive(Clone, Debug, Serialize)]
pub struct StabilityGrid {
    pub re_min: f64,
    pub re_max: f64,
    pub im_min: f64,
    pub im_max: f64,
    pub width: usize,
    pub height: usize,
    /// `|R(z)|` row major, rows running from `im_min` to `im_max`.
    pub magnitude: Vec<f64>,
}

/// Sample `|R|` on a rectangle of the complex plane.
pub fn sample_region(
    evaluate: impl Fn(Complex) -> f64,
    re: (f64, f64),
    im: (f64, f64),
    width: usize,
    height: usize,
) -> StabilityGrid {
    let mut magnitude = Vec::with_capacity(width * height);
    for row in 0..height {
        let y = im.0 + (im.1 - im.0) * row as f64 / (height - 1).max(1) as f64;
        for column in 0..width {
            let x = re.0 + (re.1 - re.0) * column as f64 / (width - 1).max(1) as f64;
            magnitude.push(evaluate(Complex::new(x, y)));
        }
    }
    StabilityGrid {
        re_min: re.0,
        re_max: re.1,
        im_min: im.0,
        im_max: im.1,
        width,
        height,
        magnitude,
    }
}

// ---------------------------------------------------------------------------
// Linear multistep methods
// ---------------------------------------------------------------------------

/// The generating polynomials of a multistep method,
/// `rho(w) = sum_j alpha_j w^(k-j)` and `sigma(w) = sum_j beta_j w^(k-j)`.
#[derive(Clone, Debug, Serialize)]
pub struct GeneratingPolynomials {
    /// Index `i` multiplies `w^i`.
    pub rho: Vec<f64>,
    pub sigma: Vec<f64>,
}

impl GeneratingPolynomials {
    pub fn from_coefficients(coefficients: &LmmCoefficients) -> GeneratingPolynomials {
        let k = coefficients.alpha.len() - 1;
        let mut rho = vec![0.0; k + 1];
        let mut sigma = vec![0.0; k + 1];
        for j in 0..=k {
            rho[k - j] = coefficients.alpha[j];
            sigma[k - j] = coefficients.beta[j];
        }
        GeneratingPolynomials { rho, sigma }
    }

    /// The boundary locus `z(theta) = rho(e^{i theta}) / sigma(e^{i theta})`.
    pub fn boundary_locus(&self, samples: usize) -> Vec<Complex> {
        (0..samples)
            .map(|i| {
                let theta = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
                let w = Complex::expi(theta);
                eval_poly(&self.rho, w) / eval_poly(&self.sigma, w)
            })
            .collect()
    }

    /// Root condition for a given `z`: every root of `rho - z sigma` inside the
    /// closed unit disc.
    pub fn is_stable_at(&self, z: Complex) -> bool {
        let degree = self.rho.len().max(self.sigma.len());
        let coefficients: Vec<Complex> = (0..degree)
            .map(|i| {
                let r = self.rho.get(i).copied().unwrap_or(0.0);
                let s = self.sigma.get(i).copied().unwrap_or(0.0);
                Complex::real(r) - z * Complex::real(s)
            })
            .collect();
        poly_roots(&coefficients)
            .into_iter()
            .all(|w| w.abs() <= 1.0 + 1e-10)
    }

    /// Zero stability, the root condition at `z = 0`.
    pub fn is_zero_stable(&self) -> bool {
        let roots = poly_roots(&self.rho.iter().map(|v| Complex::real(*v)).collect::<Vec<_>>());
        // Roots on the unit circle must be simple.
        for (i, root) in roots.iter().enumerate() {
            if root.abs() > 1.0 + 1e-8 {
                return false;
            }
            if root.abs() > 1.0 - 1e-8 {
                for (j, other) in roots.iter().enumerate() {
                    if i != j && (*root - *other).abs() < 1e-6 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Half angle of the largest wedge around the negative real axis contained
    /// in the stability region, in degrees. Ninety means A-stable.
    ///
    /// Reported to two decimals so it can be held against the published values,
    /// which for BDF3 to BDF6 are 86.03, 73.35, 51.84 and 17.84 degrees.
    pub fn alpha_angle(&self) -> f64 {
        let radii: Vec<f64> = (0..61).map(|i| 10f64.powf(-5.0 + i as f64 / 6.0)).collect();
        let wedge_is_stable = |degrees: f64| {
            let angle = std::f64::consts::PI * (1.0 - degrees / 180.0);
            radii.iter().all(|r| {
                let z = Complex::new(r * angle.cos(), r * angle.sin());
                self.is_stable_at(z) && self.is_stable_at(z.conj())
            })
        };

        if wedge_is_stable(90.0) {
            return 90.0;
        }
        let mut low = 0.0f64;
        let mut high = 90.0f64;
        if !wedge_is_stable(low) {
            return 0.0;
        }
        // Bisect on the wedge half angle. The region is a wedge by definition,
        // so stability is monotone in the angle and bisection is exact.
        for _ in 0..40 {
            let mid = 0.5 * (low + high);
            if wedge_is_stable(mid) {
                low = mid;
            } else {
                high = mid;
            }
            if high - low < 5e-3 {
                break;
            }
        }
        (low * 100.0).round() / 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn function(numerator: Vec<f64>, denominator: Vec<f64>) -> StabilityFunction {
        StabilityFunction {
            numerator,
            denominator,
        }
    }

    #[test]
    fn backward_euler_is_l_stable() {
        // R(z) = 1 / (1 - z)
        let r = function(vec![1.0], vec![1.0, -1.0]);
        assert!(r.is_a_stable());
        assert!(r.is_l_stable());
        assert_eq!(r.at_infinity(), 0.0);
        assert_eq!(r.order_of_consistency(6), 1);
    }

    #[test]
    fn trapezoidal_is_a_stable_but_not_l_stable() {
        // R(z) = (1 + z/2) / (1 - z/2)
        let r = function(vec![1.0, 0.5], vec![1.0, -0.5]);
        assert!(r.is_a_stable());
        assert!(!r.is_l_stable());
        assert!((r.at_infinity() + 1.0).abs() < 1e-12);
        assert_eq!(r.order_of_consistency(6), 2);
    }

    #[test]
    fn forward_euler_is_not_a_stable_and_reaches_minus_two() {
        let r = function(vec![1.0, 1.0], vec![1.0]);
        assert!(!r.is_a_stable());
        assert!((r.real_stability_limit() + 2.0).abs() < 1e-6);
        assert_eq!(r.order_of_consistency(6), 1);
    }

    #[test]
    fn characteristic_polynomial_of_a_known_matrix() {
        // [[1,2],[3,4]] has charpoly lambda^2 - 5 lambda - 2.
        let a = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let c = characteristic_polynomial(&a);
        assert!((c[2] - 1.0).abs() < 1e-12);
        assert!((c[1] + 5.0).abs() < 1e-12);
        assert!((c[0] + 2.0).abs() < 1e-12);
    }
}

/// Largest negative real `x` with `[x, 0]` inside the stability region.
pub fn scan_real_limit(stable: impl Fn(f64) -> bool) -> f64 {
    scan_limit(stable, -1.0)
}

impl GeneratingPolynomials {
    /// Largest root modulus of `rho - z sigma`, the multistep analogue of
    /// `|R(z)|`. Values at or below one mean the mode does not grow.
    pub fn root_radius(&self, z: Complex) -> f64 {
        let degree = self.rho.len().max(self.sigma.len());
        let coefficients: Vec<Complex> = (0..degree)
            .map(|i| {
                let r = self.rho.get(i).copied().unwrap_or(0.0);
                let s = self.sigma.get(i).copied().unwrap_or(0.0);
                Complex::real(r) - z * Complex::real(s)
            })
            .collect();
        poly_roots(&coefficients)
            .into_iter()
            .fold(0.0f64, |acc, w| acc.max(w.abs()))
    }

    /// The part of the boundary locus that actually bounds the stability
    /// region, with everything else replaced by `NaN`.
    ///
    /// Every point of the locus has a root on the unit circle, which is why the
    /// curve is where the region can end. It is not where it does end: on the
    /// outer loops of a high order Adams or BDF family another root is already
    /// outside the circle, so the curve runs through the unstable set without
    /// bounding anything. Keeping only the points whose largest root is the one
    /// on the circle leaves exactly the boundary, and does it in closed form
    /// where a contour through a sampled grid would have to guess.
    pub fn region_boundary(&self, samples: usize) -> Vec<Complex> {
        self.boundary_locus(samples)
            .into_iter()
            .map(|z| {
                if z.abs().is_finite() && self.root_radius(z) <= 1.0 + 1e-6 {
                    z
                } else {
                    Complex::new(f64::NAN, f64::NAN)
                }
            })
            .collect()
    }
}
