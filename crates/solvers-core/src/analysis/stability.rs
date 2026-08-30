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
//! infinity and a test for A-stability on the imaginary axis, instead of a
//! picture that has to be eyeballed.
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
use crate::method::{LmmCoefficients, RkTableau, RosenbrockTableau};
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
        let a = tableau.a.map(|v| v.value());
        let b: Vec<f64> = tableau.b.iter().map(|v| v.value()).collect();
        StabilityFunction::from_parts(&a, &b)
    }

    /// Build `R(z)` for a Rosenbrock method.
    ///
    /// On a linear problem the Jacobian is the problem, so the two couplings
    /// add and the method behaves exactly like the diagonally implicit
    /// Runge-Kutta method with `A = alpha + gamma`. That is the reason a ROW
    /// method has a DIRK stability function and can be L-stable at all.
    pub fn from_rosenbrock(tableau: &RosenbrockTableau) -> StabilityFunction {
        let a = tableau.implicit_matrix().map(|v| v.value());
        let b: Vec<f64> = tableau.b.iter().map(|v| v.value()).collect();
        StabilityFunction::from_parts(&a, &b)
    }

    /// `R(z) = det(I - zA + z e b^T) / det(I - zA)` from the two pieces it is
    /// made of.
    pub fn from_parts(a: &Matrix<f64>, b: &[f64]) -> StabilityFunction {
        let s = b.len();

        // det(I - zA) = z^s * charpoly_A(1/z), so the coefficients just reverse.
        let ca = characteristic_polynomial(a);
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

    /// Dispersion and dissipation order.
    ///
    /// On the imaginary axis the exact solution turns without growing, and a
    /// method does neither exactly. Writing `W(y) = R(iy) exp(-iy)`, the real
    /// part of its first nonzero term is the amplitude it gains or loses per
    /// step and the imaginary part is the phase it runs early or late, so the
    /// two orders are read straight off one series.
    ///
    /// `None` for the dissipation order means there is none at any order, which
    /// is what `|R(iy)| = 1` looks like and is exactly the property that makes
    /// the Gauss methods worth their cost on an oscillator.
    ///
    /// Reference: P. J. van der Houwen, B. P. Sommeijer, "Explicit Runge-Kutta
    /// (-Nystroem) methods with reduced phase errors for computing oscillating
    /// solutions", SIAM J. Numer. Anal. 24(3), 1987, doi:10.1137/0724041
    pub fn phase_and_amplitude_order(&self, max_order: usize) -> (Option<usize>, Option<usize>) {
        let terms = max_order + 2;
        let taylor = self.taylor(terms);

        let mut factorial = vec![1.0f64; terms + 1];
        for k in 1..=terms {
            factorial[k] = factorial[k - 1] * k as f64;
        }
        let power_of_i = |k: usize| match k % 4 {
            0 => Complex::new(1.0, 0.0),
            1 => Complex::new(0.0, 1.0),
            2 => Complex::new(-1.0, 0.0),
            _ => Complex::new(0.0, -1.0),
        };

        // `W = R(iy) exp(-iy)`, term by term.
        let mut w = vec![Complex::new(0.0, 0.0); terms + 1];
        for (m, slot) in w.iter_mut().enumerate() {
            for k in 0..=m {
                let left = power_of_i(k) * Complex::real(taylor[k]);
                let right = power_of_i(m - k).conj() / Complex::real(factorial[m - k]);
                *slot = *slot + left * right;
            }
        }

        // A term of `W` is a sum of about `m` products of size `2^m / m!`, so
        // that is the size of the round off in it. Anything a few orders above
        // is a real coefficient, anything below is the arithmetic.
        let significant = |m: usize| 1e-7 * 2f64.powi(m as i32) / factorial[m];

        // The phase is the argument, whose first term is the imaginary part.
        let dispersion = (1..=terms)
            .find(|m| w[*m].im.abs() > significant(*m))
            .map(|m| m - 1);

        // The amplitude is the modulus, and the modulus is not the real part.
        // A method that only runs early or late has `W` on the unit circle with
        // a real part that departs from one at twice the order of its phase
        // error, and reading that as an amplitude error would credit a Gauss
        // method, whose amplitude is exactly one, with a dissipation of order
        // thirteen. So `|W|^2` is formed and its first term taken.
        let mut squared = vec![0.0f64; terms + 1];
        for (n, slot) in squared.iter_mut().enumerate() {
            for k in 0..=n {
                let product = w[k] * w[n - k].conj();
                *slot += product.re;
            }
        }
        let dissipation = (1..=terms)
            .find(|n| squared[*n].abs() > significant(*n))
            .map(|n| n - 1);

        (dispersion, dissipation)
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
    ///
    /// The condition on the axis is `|N(iy)| <= |D(iy)|`, and it is decided by
    /// evaluating the two rather than by expanding their difference. That
    /// expansion is the E-polynomial, and from six stages upwards its
    /// coefficients are differences of numbers that agree to fifteen digits:
    /// for Radau IIA, where the truth is a single positive term in `u^s`, what
    /// comes out is round off whose leading coefficients have arbitrary signs.
    /// The roots of it are still worth having, because they say where `|R|`
    /// could touch one and therefore where a crossing would hide.
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

        let mut probes: Vec<f64> = (0..=240).map(|i| 10f64.powf(-8.0 + i as f64 / 15.0)).collect();
        probes.push(0.0);
        for square in positive_real_roots(&self.e_polynomial()) {
            let y = square.sqrt();
            probes.extend([0.98 * y, y, 1.02 * y]);
        }
        probes.into_iter().all(|y| {
            let modulus = self.eval(Complex::new(0.0, y)).abs();
            modulus.is_finite() && modulus <= 1.0 + 1e-9
        })
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

/// The positive real roots of a real polynomial, in increasing order.
fn positive_real_roots(p: &[f64]) -> Vec<f64> {
    let mut trimmed = p.to_vec();
    while trimmed.len() > 1 && trimmed.last().map_or(false, |v| v.abs() < 1e-13) {
        trimmed.pop();
    }
    if trimmed.len() < 2 {
        return Vec::new();
    }
    let coefficients: Vec<Complex> = trimmed.iter().map(|v| Complex::real(*v)).collect();
    let mut roots: Vec<f64> = poly_roots(&coefficients)
        .into_iter()
        .filter(|r| r.im.abs() < 1e-7 && r.re > 0.0)
        .map(|r| r.re)
        .collect();
    roots.sort_by(f64::total_cmp);
    roots
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

/// Trace the curve where `value` crosses `0`, by adaptive subdivision.
///
/// A uniform grid spends almost everything it has where nothing happens. The
/// boundary of a stability region is a curve, so the samples belong on the
/// curve: a coarse grid is laid down, and only the cells whose corners disagree
/// about which side of the level they are on get subdivided, recursively. The
/// cost then grows with the length of the boundary rather than with the area of
/// the window, which is what makes an effective resolution of a thousand points
/// across affordable at all.
///
/// The caller passes a logarithm rather than the modulus itself. It is finite
/// at a pole, it is far closer to linear across a cell, and a crossing found by
/// interpolating it lands where the crossing actually is.
///
/// Returns the segments end to end, each pair separated by a `NaN` point so a
/// plot breaks the line between them.
pub fn trace_zero_level(
    value: impl Fn(Complex) -> f64,
    re: (f64, f64),
    im: (f64, f64),
    coarse: usize,
    depth: u32,
) -> Vec<Complex> {
    let mut out = Vec::new();
    let mut segments = Vec::new();
    let coarse = coarse.max(2);
    let dx = (re.1 - re.0) / coarse as f64;
    let dy = (im.1 - im.0) / coarse as f64;

    // One row of corner values at a time, so the coarse grid is evaluated once
    // rather than four times over.
    let corner = |i: usize, j: usize| Complex::new(re.0 + i as f64 * dx, im.0 + j as f64 * dy);
    let mut lower: Vec<f64> = (0..=coarse).map(|i| value(corner(i, 0))).collect();
    for j in 0..coarse {
        let upper: Vec<f64> = (0..=coarse).map(|i| value(corner(i, j + 1))).collect();
        for i in 0..coarse {
            let cell = Cell {
                x0: re.0 + i as f64 * dx,
                y0: im.0 + j as f64 * dy,
                x1: re.0 + (i + 1) as f64 * dx,
                y1: im.0 + (j + 1) as f64 * dy,
                v00: lower[i],
                v10: lower[i + 1],
                v01: upper[i],
                v11: upper[i + 1],
            };
            subdivide(&value, cell, depth, &mut segments);
        }
        lower = upper;
    }

    for chain in link(&segments) {
        out.extend(chain);
        out.push(Complex::new(f64::NAN, f64::NAN));
    }
    out
}

/// Chain the segments into polylines.
///
/// Marching squares produces the curve as unordered pieces, and a plot handed
/// those draws one path element per piece: a thousand nodes in the document for
/// one curve. Cells at the same depth compute their shared edge crossing from
/// the same two corner values, so the pieces meet exactly and can be chained by
/// looking an endpoint up rather than by searching for a near match.
fn link(segments: &[(Complex, Complex)]) -> Vec<Vec<Complex>> {
    use std::collections::HashMap;

    // Positive and negative zero are the same point and have to hash alike.
    let key = |z: &Complex| (((z.re + 0.0).to_bits()), ((z.im + 0.0).to_bits()));

    let mut at: HashMap<(u64, u64), Vec<usize>> = HashMap::new();
    for (index, (a, b)) in segments.iter().enumerate() {
        at.entry(key(a)).or_default().push(index);
        at.entry(key(b)).or_default().push(index);
    }

    let mut used = vec![false; segments.len()];
    let mut chains = Vec::new();
    // Follow the curve from one end of a piece, taking whichever unused piece
    // starts where the last one stopped.
    let step = |end: &Complex, used: &mut Vec<bool>| -> Option<Complex> {
        for &candidate in at.get(&key(end))? {
            if used[candidate] {
                continue;
            }
            let (a, b) = segments[candidate];
            let next = if key(&a) == key(end) { b } else { a };
            used[candidate] = true;
            return Some(next);
        }
        None
    };

    for start in 0..segments.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let (a, b) = segments[start];
        let mut chain = vec![a, b];
        while let Some(next) = step(&chain[chain.len() - 1], &mut used) {
            chain.push(next);
        }
        // And backwards from the other end, since the piece this started from
        // need not have been the first one on the curve. Collected in reverse
        // and turned around once, rather than pushed onto the front each time.
        let mut head = Vec::new();
        let mut front = chain[0];
        while let Some(previous) = step(&front, &mut used) {
            head.push(previous);
            front = previous;
        }
        head.reverse();
        head.extend(chain);
        chains.push(head);
    }
    chains
}

#[derive(Copy, Clone)]
struct Cell {
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    /// Corner values, named by their position: `v01` is left and top.
    v00: f64,
    v10: f64,
    v01: f64,
    v11: f64,
}

impl Cell {
    /// Whether the four corners agree about which side of zero they are on.
    fn uniform(&self) -> bool {
        let positive = |v: f64| v > 0.0;
        let first = positive(self.v00);
        [self.v10, self.v01, self.v11].iter().all(|v| positive(*v) == first)
    }
}

fn subdivide(
    value: &impl Fn(Complex) -> f64,
    cell: Cell,
    depth: u32,
    out: &mut Vec<(Complex, Complex)>,
) {
    if cell.uniform() {
        return;
    }
    if depth == 0 {
        emit(cell, out);
        return;
    }
    let xm = 0.5 * (cell.x0 + cell.x1);
    let ym = 0.5 * (cell.y0 + cell.y1);
    let at = |x: f64, y: f64| value(Complex::new(x, y));
    let vs = at(xm, cell.y0);
    let vn = at(xm, cell.y1);
    let vw = at(cell.x0, ym);
    let ve = at(cell.x1, ym);
    let vc = at(xm, ym);

    let quarters = [
        Cell { x0: cell.x0, y0: cell.y0, x1: xm, y1: ym, v00: cell.v00, v10: vs, v01: vw, v11: vc },
        Cell { x0: xm, y0: cell.y0, x1: cell.x1, y1: ym, v00: vs, v10: cell.v10, v01: vc, v11: ve },
        Cell { x0: cell.x0, y0: ym, x1: xm, y1: cell.y1, v00: vw, v10: vc, v01: cell.v01, v11: vn },
        Cell { x0: xm, y0: ym, x1: cell.x1, y1: cell.y1, v00: vc, v10: ve, v01: vn, v11: cell.v11 },
    ];
    for quarter in quarters {
        subdivide(value, quarter, depth - 1, out);
    }
}

/// Marching squares on one cell, at the finest level.
fn emit(cell: Cell, out: &mut Vec<(Complex, Complex)>) {
    // Where the level crosses an edge, by linear interpolation between its two
    // ends. The corners are known to disagree or this cell would not be here.
    let cross = |a: f64, b: f64| a / (a - b);
    let mut points = Vec::with_capacity(4);
    if (cell.v00 > 0.0) != (cell.v10 > 0.0) {
        let t = cross(cell.v00, cell.v10);
        points.push(Complex::new(cell.x0 + t * (cell.x1 - cell.x0), cell.y0));
    }
    if (cell.v10 > 0.0) != (cell.v11 > 0.0) {
        let t = cross(cell.v10, cell.v11);
        points.push(Complex::new(cell.x1, cell.y0 + t * (cell.y1 - cell.y0)));
    }
    if (cell.v01 > 0.0) != (cell.v11 > 0.0) {
        let t = cross(cell.v01, cell.v11);
        points.push(Complex::new(cell.x0 + t * (cell.x1 - cell.x0), cell.y1));
    }
    if (cell.v00 > 0.0) != (cell.v01 > 0.0) {
        let t = cross(cell.v00, cell.v01);
        points.push(Complex::new(cell.x0, cell.y0 + t * (cell.y1 - cell.y0)));
    }
    // Two crossings make one segment. Four make a saddle, and at this size
    // either pairing is within a pixel of the other.
    let mut index = 0;
    while index + 1 < points.len() {
        out.push((points[index], points[index + 1]));
        index += 2;
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

    /// Explicit Euler has `R(z) = 1 + z`, so its boundary is the circle of
    /// radius one about minus one, exactly. Every point the tracer returns has
    /// to sit on it, which is the only check that says the adaptive
    /// subdivision, the interpolation and the cell wiring are all right.
    #[test]
    fn the_traced_boundary_lands_on_the_curve_it_is_tracing() {
        let function = StabilityFunction {
            numerator: vec![1.0, 1.0],
            denominator: vec![1.0],
        };
        let curve = trace_zero_level(
            |z| {
                let magnitude = function.eval(z).abs();
                if magnitude > 0.0 {
                    magnitude.log10().clamp(-30.0, 30.0)
                } else {
                    -30.0
                }
            },
            (-3.0, 1.0),
            (-2.0, 2.0),
            56,
            4,
        );

        let points: Vec<&Complex> = curve.iter().filter(|z| z.re.is_finite()).collect();
        assert!(points.len() > 200, "only {} points on the boundary", points.len());
        let worst = points
            .iter()
            .map(|z| ((z.re + 1.0).hypot(z.im) - 1.0).abs())
            .fold(0.0f64, f64::max);
        assert!(worst < 2e-4, "worst departure from the unit circle {worst:.2e}");
    }


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
