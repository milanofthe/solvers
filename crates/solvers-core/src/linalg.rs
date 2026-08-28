//! Dense linear algebra, generic over the scalar field.
//!
//! The same LU is used three times: to solve the Newton systems during time
//! stepping (`f64`), to evaluate stability functions on the complex plane
//! (`Complex`) and to determine variable step multistep coefficients.

use crate::num::Field;

#[derive(Clone, Debug, PartialEq)]
pub struct Matrix<T> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
}

impl<T: Field> Matrix<T> {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Matrix {
            rows,
            cols,
            data: vec![T::zero(); rows * cols],
        }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = Matrix::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = T::one();
        }
        m
    }

    pub fn from_rows(rows: &[Vec<T>]) -> Self {
        let r = rows.len();
        let c = rows.first().map_or(0, |row| row.len());
        let mut m = Matrix::zeros(r, c);
        for (i, row) in rows.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                m[(i, j)] = *v;
            }
        }
        m
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn row(&self, i: usize) -> &[T] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }

    pub fn fill(&mut self, v: T) {
        for x in self.data.iter_mut() {
            *x = v;
        }
    }

    /// Map every entry into another field, e.g. rational tableau to complex.
    pub fn map<U: Field>(&self, f: impl Fn(T) -> U) -> Matrix<U> {
        Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().map(|v| f(*v)).collect(),
        }
    }

    /// Matrix vector product.
    pub fn mul_vec(&self, x: &[T], out: &mut [T]) {
        for i in 0..self.rows {
            let mut acc = T::zero();
            for j in 0..self.cols {
                acc = acc + self[(i, j)] * x[j];
            }
            out[i] = acc;
        }
    }
}

impl<T> std::ops::Index<(usize, usize)> for Matrix<T> {
    type Output = T;
    fn index(&self, (i, j): (usize, usize)) -> &T {
        &self.data[i * self.cols + j]
    }
}

impl<T> std::ops::IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        &mut self.data[i * self.cols + j]
    }
}

/// LU factorization with partial pivoting.
#[derive(Clone, Debug)]
pub struct Lu<T> {
    lu: Matrix<T>,
    perm: Vec<usize>,
    sign: f64,
    singular: bool,
}

impl<T: Field> Lu<T> {
    /// Factor a square matrix in place.
    ///
    /// A singular pivot is recorded rather than rejected so that callers doing
    /// determinant style work (stability functions at poles) can still ask.
    pub fn factor(mut a: Matrix<T>) -> Lu<T> {
        let n = a.rows();
        debug_assert_eq!(n, a.cols(), "LU requires a square matrix");
        let mut perm: Vec<usize> = (0..n).collect();
        let mut sign = 1.0;
        let mut singular = false;

        for k in 0..n {
            let mut pivot = k;
            let mut best = a[(k, k)].magnitude();
            for i in (k + 1)..n {
                let m = a[(i, k)].magnitude();
                if m > best {
                    best = m;
                    pivot = i;
                }
            }
            if best == 0.0 || !best.is_finite() {
                singular = true;
                continue;
            }
            if pivot != k {
                for j in 0..n {
                    let tmp = a[(k, j)];
                    a[(k, j)] = a[(pivot, j)];
                    a[(pivot, j)] = tmp;
                }
                perm.swap(k, pivot);
                sign = -sign;
            }
            let d = a[(k, k)];
            for i in (k + 1)..n {
                let factor = a[(i, k)] / d;
                a[(i, k)] = factor;
                if factor.is_zero() {
                    continue;
                }
                for j in (k + 1)..n {
                    let v = a[(k, j)];
                    a[(i, j)] = a[(i, j)] - factor * v;
                }
            }
        }

        Lu {
            lu: a,
            perm,
            sign,
            singular,
        }
    }

    pub fn is_singular(&self) -> bool {
        self.singular
    }

    pub fn dim(&self) -> usize {
        self.lu.rows()
    }

    /// Solve `A x = b` in place, `b` is overwritten with the solution.
    pub fn solve_in_place(&self, b: &mut [T]) -> bool {
        if self.singular {
            return false;
        }
        let n = self.dim();
        let mut y = vec![T::zero(); n];
        for i in 0..n {
            y[i] = b[self.perm[i]];
        }
        // Forward substitution, unit lower triangular.
        for i in 1..n {
            let mut acc = y[i];
            for j in 0..i {
                acc = acc - self.lu[(i, j)] * y[j];
            }
            y[i] = acc;
        }
        // Back substitution.
        for i in (0..n).rev() {
            let mut acc = y[i];
            for j in (i + 1)..n {
                acc = acc - self.lu[(i, j)] * y[j];
            }
            y[i] = acc / self.lu[(i, i)];
        }
        b.copy_from_slice(&y);
        true
    }

    pub fn solve(&self, b: &[T]) -> Option<Vec<T>> {
        let mut x = b.to_vec();
        if self.solve_in_place(&mut x) {
            Some(x)
        } else {
            None
        }
    }

    pub fn determinant(&self) -> T {
        if self.singular {
            return T::zero();
        }
        let mut det = T::from_f64(self.sign);
        for i in 0..self.dim() {
            det = det * self.lu[(i, i)];
        }
        det
    }
}

/// Solve a small dense system directly.
pub fn solve(a: Matrix<f64>, b: &[f64]) -> Option<Vec<f64>> {
    Lu::factor(a).solve(b)
}

/// Roots of a polynomial by Durand-Kerner iteration.
///
/// `coeffs[i]` multiplies `x^i`. Used for the root condition of linear
/// multistep methods, where the degree is small and the accuracy needed is
/// only enough to compare moduli against one.
pub fn poly_roots(coeffs: &[crate::num::Complex]) -> Vec<crate::num::Complex> {
    use crate::num::Complex;

    let mut c: Vec<Complex> = coeffs.to_vec();
    while c.len() > 1 && c.last().map_or(false, |v| v.abs() < 1e-14) {
        c.pop();
    }
    let degree = c.len().saturating_sub(1);
    if degree == 0 {
        return Vec::new();
    }
    let lead = c[degree];
    let monic: Vec<Complex> = c.iter().map(|v| *v / lead).collect();

    // Spread the initial guesses on a circle to avoid symmetric stalling.
    let seed = Complex::new(0.4, 0.9);
    let mut roots: Vec<Complex> = (0..degree)
        .map(|i| {
            let mut p = Complex::new(1.0, 0.0);
            for _ in 0..i {
                p = p * seed;
            }
            p
        })
        .collect();

    let eval = |x: Complex| -> Complex {
        let mut acc = Complex::new(0.0, 0.0);
        for coeff in monic.iter().rev() {
            acc = acc * x + *coeff;
        }
        acc
    };

    for _ in 0..200 {
        let mut delta = 0.0f64;
        for i in 0..degree {
            let mut denom = Complex::new(1.0, 0.0);
            for j in 0..degree {
                if i != j {
                    denom = denom * (roots[i] - roots[j]);
                }
            }
            if denom.abs() < 1e-300 {
                continue;
            }
            let step = eval(roots[i]) / denom;
            roots[i] = roots[i] - step;
            delta = delta.max(step.abs());
        }
        if delta < 1e-14 {
            break;
        }
    }
    roots
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::num::Complex;

    #[test]
    fn lu_solves_a_real_system() {
        let a = Matrix::from_rows(&[vec![2.0, 1.0, -1.0], vec![-3.0, -1.0, 2.0], vec![-2.0, 1.0, 2.0]]);
        let x = solve(a, &[8.0, -11.0, -3.0]).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-12);
        assert!((x[1] - 3.0).abs() < 1e-12);
        assert!((x[2] + 1.0).abs() < 1e-12);
    }

    #[test]
    fn determinant_matches() {
        let a = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert!((Lu::factor(a).determinant() + 2.0).abs() < 1e-12);
    }

    #[test]
    fn complex_system() {
        let a = Matrix::from_rows(&[
            vec![Complex::new(1.0, 1.0), Complex::new(0.0, 0.0)],
            vec![Complex::new(0.0, 0.0), Complex::new(2.0, 0.0)],
        ]);
        let x = Lu::factor(a).solve(&[Complex::new(2.0, 0.0), Complex::new(4.0, 0.0)]).unwrap();
        assert!((x[0] - Complex::new(1.0, -1.0)).abs() < 1e-12);
        assert!((x[1] - Complex::new(2.0, 0.0)).abs() < 1e-12);
    }

    #[test]
    fn roots_of_a_quadratic() {
        // x^2 - 3x + 2 has roots 1 and 2.
        let roots = poly_roots(&[
            Complex::real(2.0),
            Complex::real(-3.0),
            Complex::real(1.0),
        ]);
        let mut moduli: Vec<f64> = roots.iter().map(|r| r.re).collect();
        moduli.sort_by(|a, b| a.total_cmp(b));
        assert!((moduli[0] - 1.0).abs() < 1e-8);
        assert!((moduli[1] - 2.0).abs() < 1e-8);
    }
}
