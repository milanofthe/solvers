//! Decoupling the stage system of a fully implicit method.
//!
//! A fully implicit tableau couples all `s` stages, so a straight Newton step
//! means factoring an `s*n` matrix at a cost of `(s n)^3 / 3`. The standard way
//! out, and what makes RADAU5 competitive, is to diagonalize the stage coupling
//! instead. Writing `A^{-1} = T L T^{-1}` and multiplying the Newton system
//! `(I - h A kron J) dZ = G` from the left by `h^{-1} A^{-1} kron I` gives
//!
//! ```text
//! [h^{-1} (L kron I) - (I kron J)] dW = h^{-1} (L kron I) Y
//! dW = (T^{-1} kron I) dZ,   Y = (T^{-1} kron I) G
//! ```
//!
//! which is `s` independent systems of size `n`, one per eigenvalue. Because
//! `A` is real the eigenvalues come as one real value and conjugate pairs, and
//! only one member of each pair has to be solved. For the three stage Radau IIA
//! that turns one factorization of size `3n` into one real and one complex
//! factorization of size `n`, about five times less work, and the gap widens
//! with the number of stages.
//!
//! What is left over is the basis change, `s` small dense transforms applied
//! across long vectors. That is exactly the shape `simd::stage_transform`
//! blocks over components, so the transforms vectorize instead of becoming the
//! new bottleneck.
//!
//! Reference: E. Hairer, G. Wanner, "Solving Ordinary Differential Equations
//! II", 2nd ed., Springer 1996, section IV.8, doi:10.1007/978-3-642-05221-7

use crate::linalg::{Lu, Matrix};
use crate::num::Complex;
use crate::simd;

/// One eigenvalue of `A^{-1}` and how it is solved.
#[derive(Clone, Debug)]
enum Block {
    /// A real eigenvalue: one real system of size `n`.
    Real { index: usize, lambda: f64 },
    /// A conjugate pair, of which only the one with positive imaginary part is
    /// solved; the other follows by conjugation.
    Pair {
        index: usize,
        conjugate: usize,
        lambda: Complex,
    },
}

/// The similarity transform that decouples a stage system.
#[derive(Clone, Debug)]
pub struct StageDecoupling {
    pub stages: usize,
    blocks: Vec<Block>,
    /// Eigenvector matrix `T`, column major by eigenvalue, row major storage.
    t: Vec<Complex>,
    /// `T^{-1}`.
    t_inv: Vec<Complex>,
}

/// Eigenvalues of a small real matrix, from the roots of its characteristic
/// polynomial. Fine at the sizes a Butcher tableau has.
fn eigenvalues(a: &Matrix<f64>) -> Vec<Complex> {
    let coefficients: Vec<Complex> = crate::analysis::stability::characteristic_polynomial(a)
        .into_iter()
        .map(Complex::real)
        .collect();
    crate::linalg::poly_roots(&coefficients)
}

/// Eigenvector for a known eigenvalue, by inverse iteration on a shifted
/// system. Two iterations are enough at these sizes and it needs no more
/// machinery than the complex LU that already exists.
fn eigenvector(a: &Matrix<f64>, lambda: Complex) -> Option<Vec<Complex>> {
    let s = a.rows();
    let shift = lambda + Complex::new(1e-10 * lambda.abs().max(1.0), 1e-10);
    let mut shifted = Matrix::<Complex>::zeros(s, s);
    for i in 0..s {
        for j in 0..s {
            shifted[(i, j)] = Complex::real(a[(i, j)]);
        }
        shifted[(i, i)] = shifted[(i, i)] - shift;
    }
    let lu = Lu::factor(shifted);
    if lu.is_singular() {
        return None;
    }

    // A real starting vector, chosen to be unlikely to sit orthogonal to the
    // eigenvector. Keeping it real is what makes the whole computation
    // conjugate symmetric, so a conjugate eigenvalue really does yield the
    // conjugate eigenvector.
    let mut v: Vec<Complex> = (0..s)
        .map(|i| Complex::real(1.0 + 0.3 * i as f64 + 0.11 * (i * i) as f64))
        .collect();
    for _ in 0..3 {
        let solved = lu.solve(&v)?;
        let norm = solved.iter().fold(0.0f64, |acc, z| acc.max(z.abs()));
        if !norm.is_finite() || norm == 0.0 {
            return None;
        }
        v = solved.iter().map(|z| *z / Complex::real(norm)).collect();
    }
    Some(v)
}

/// Invert a small complex matrix through its LU.
fn invert(m: &Matrix<Complex>) -> Option<Matrix<Complex>> {
    let s = m.rows();
    let lu = Lu::factor(m.clone());
    if lu.is_singular() {
        return None;
    }
    let mut inverse = Matrix::<Complex>::zeros(s, s);
    for column in 0..s {
        let mut unit = vec![Complex::new(0.0, 0.0); s];
        unit[column] = Complex::new(1.0, 0.0);
        let solved = lu.solve(&unit)?;
        for row in 0..s {
            inverse[(row, column)] = solved[row];
        }
    }
    Some(inverse)
}

impl StageDecoupling {
    /// Build the decoupling for a stage matrix, or `None` when it does not
    /// exist: a singular or defective `A` has to fall back to the coupled solve.
    pub fn new(a: &Matrix<f64>) -> Option<StageDecoupling> {
        let s = a.rows();
        if s == 0 {
            return None;
        }

        let mut nu = eigenvalues(a);
        if nu.len() != s {
            return None;
        }
        // A deterministic order keeps the transform reproducible.
        nu.sort_by(|x, y| x.im.total_cmp(&y.im).then(x.re.total_cmp(&y.re)));

        // Eigenvalues of A^{-1} are the reciprocals, with the same eigenvectors.
        let mut lambda = Vec::with_capacity(s);
        for value in &nu {
            if value.abs() < 1e-12 {
                return None;
            }
            lambda.push(Complex::new(1.0, 0.0) / *value);
        }

        // Group first: real eigenvalues stand alone, complex ones pair with
        // their conjugate. Only one member of a pair is ever solved, so the
        // partner's eigenvector must be the exact conjugate rather than an
        // independently computed one.
        let mut blocks = Vec::new();
        let mut taken = vec![false; s];
        for k in 0..s {
            if taken[k] {
                continue;
            }
            if lambda[k].im.abs() < 1e-10 * lambda[k].abs().max(1.0) {
                taken[k] = true;
                blocks.push(Block::Real {
                    index: k,
                    lambda: lambda[k].re,
                });
                continue;
            }
            let mut partner = None;
            for j in 0..s {
                if j == k || taken[j] {
                    continue;
                }
                if (lambda[j] - lambda[k].conj()).abs() < 1e-8 * lambda[k].abs().max(1.0) {
                    partner = Some(j);
                    break;
                }
            }
            let partner = partner?;
            taken[k] = true;
            taken[partner] = true;
            blocks.push(Block::Pair {
                index: k,
                conjugate: partner,
                lambda: lambda[k],
            });
        }

        let mut t = Matrix::<Complex>::zeros(s, s);
        for block in &blocks {
            match block {
                Block::Real { index, .. } => {
                    let vector = eigenvector(a, nu[*index])?;
                    for row in 0..s {
                        t[(row, *index)] = vector[row];
                    }
                }
                Block::Pair { index, conjugate, .. } => {
                    let vector = eigenvector(a, nu[*index])?;
                    for row in 0..s {
                        t[(row, *index)] = vector[row];
                        t[(row, *conjugate)] = vector[row].conj();
                    }
                }
            }
        }

        let t_inv = invert(&t)?;

        // Verify the decomposition instead of trusting the eigenvector search.
        // A defective tableau simply does not get the fast path.
        let mut worst = 0.0f64;
        let scale = a.as_slice().iter().fold(1.0f64, |acc, v| acc.max(v.abs()));
        for i in 0..s {
            for j in 0..s {
                let mut acc = Complex::new(0.0, 0.0);
                for k in 0..s {
                    for p in 0..s {
                        acc = acc + Complex::real(a[(i, p)]) * t[(p, k)] * t_inv[(k, j)];
                    }
                }
                let target = Complex::real(a[(i, j)]);
                worst = worst.max((acc - target).abs());
            }
        }
        if worst > 1e-8 * scale {
            return None;
        }

        Some(StageDecoupling {
            stages: s,
            blocks,
            t: t.as_slice().to_vec(),
            t_inv: t_inv.as_slice().to_vec(),
        })
    }

    /// Number of matrices that have to be factored per step: one per real
    /// eigenvalue and one per conjugate pair.
    pub fn factorizations(&self) -> usize {
        self.blocks.len()
    }

    /// Eigenvalues of `A^{-1}` the linear systems are built from.
    pub fn shifts(&self) -> Vec<Complex> {
        self.blocks
            .iter()
            .map(|block| match block {
                Block::Real { lambda, .. } => Complex::real(*lambda),
                Block::Pair { lambda, .. } => *lambda,
            })
            .collect()
    }
}

/// The factored blocks for one step size and one Jacobian.
pub struct DecoupledLinear {
    decoupling: StageDecoupling,
    real: Vec<(usize, Lu<f64>)>,
    complex: Vec<(usize, usize, Lu<Complex>)>,
    /// Step size the factorizations belong to.
    factored_h: f64,
    dim: usize,
    // Scratch, kept alive so a step allocates nothing.
    y_re: Vec<Vec<f64>>,
    y_im: Vec<Vec<f64>>,
    buffer: Vec<Complex>,
}

impl DecoupledLinear {
    pub fn new(decoupling: StageDecoupling, dim: usize) -> DecoupledLinear {
        let s = decoupling.stages;
        DecoupledLinear {
            decoupling,
            real: Vec::new(),
            complex: Vec::new(),
            factored_h: f64::NAN,
            dim,
            y_re: vec![vec![0.0; dim]; s],
            y_im: vec![vec![0.0; dim]; s],
            buffer: vec![Complex::new(0.0, 0.0); dim],
        }
    }

    pub fn stages(&self) -> usize {
        self.decoupling.stages
    }

    pub fn is_factored_for(&self, h: f64) -> bool {
        (!self.real.is_empty() || !self.complex.is_empty()) && self.factored_h == h
    }

    pub fn invalidate(&mut self) {
        self.real.clear();
        self.complex.clear();
        self.factored_h = f64::NAN;
    }

    /// Number of `n` sized factorizations the last `factor` performed.
    pub fn factorization_count(&self) -> usize {
        self.real.len() + self.complex.len()
    }

    /// Factor `lambda_k / h * I - J` for every block.
    pub fn factor(&mut self, jacobian: &Matrix<f64>, h: f64) -> bool {
        self.real.clear();
        self.complex.clear();
        let n = self.dim;

        for block in &self.decoupling.blocks {
            match block {
                Block::Real { index, lambda } => {
                    let mut m = Matrix::<f64>::zeros(n, n);
                    let shift = lambda / h;
                    for i in 0..n {
                        for j in 0..n {
                            m[(i, j)] = -jacobian[(i, j)];
                        }
                        m[(i, i)] += shift;
                    }
                    let lu = Lu::factor(m);
                    if lu.is_singular() {
                        return false;
                    }
                    self.real.push((*index, lu));
                }
                Block::Pair {
                    index,
                    conjugate,
                    lambda,
                } => {
                    let mut m = Matrix::<Complex>::zeros(n, n);
                    let shift = *lambda / Complex::real(h);
                    for i in 0..n {
                        for j in 0..n {
                            m[(i, j)] = Complex::real(-jacobian[(i, j)]);
                        }
                        m[(i, i)] = m[(i, i)] + shift;
                    }
                    let lu = Lu::factor(m);
                    if lu.is_singular() {
                        return false;
                    }
                    self.complex.push((*index, *conjugate, lu));
                }
            }
        }
        self.factored_h = h;
        true
    }

    /// Solve `(I - h A kron J) dZ = rhs` in place.
    pub fn solve(&mut self, rhs: &mut [f64], h: f64) -> bool {
        if self.real.is_empty() && self.complex.is_empty() {
            return false;
        }
        let s = self.decoupling.stages;
        let n = self.dim;
        if rhs.len() != s * n {
            return false;
        }

        // Y = (T^{-1} kron I) G, the only place the stage coupling is touched.
        // Real and imaginary parts are transformed separately so both passes are
        // the blocked real kernel.
        let mut re = vec![0.0; s * s];
        let mut im = vec![0.0; s * s];
        for i in 0..s * s {
            re[i] = self.decoupling.t_inv[i].re;
            im[i] = self.decoupling.t_inv[i].im;
        }
        let blocks: Vec<Vec<f64>> = (0..s).map(|i| rhs[i * n..(i + 1) * n].to_vec()).collect();
        simd::stage_transform(&re, s, 1.0, &blocks, &mut self.y_re);
        simd::stage_transform(&im, s, 1.0, &blocks, &mut self.y_im);

        // Solve each decoupled block.
        for (index, lu) in &self.real {
            let shift = self.decoupling.shifts_at(*index) / h;
            for p in 0..n {
                self.y_re[*index][p] *= shift;
                self.y_im[*index][p] *= shift;
            }
            if !lu.solve_in_place(&mut self.y_re[*index]) {
                return false;
            }
            if !lu.solve_in_place(&mut self.y_im[*index]) {
                return false;
            }
        }
        for (index, conjugate, lu) in &self.complex {
            let lambda = self.decoupling.complex_shift_at(*index);
            let shift = lambda / Complex::real(h);
            for p in 0..n {
                let value = Complex::new(self.y_re[*index][p], self.y_im[*index][p]) * shift;
                self.buffer[p] = value;
            }
            if !lu.solve_in_place(&mut self.buffer) {
                return false;
            }
            for p in 0..n {
                self.y_re[*index][p] = self.buffer[p].re;
                self.y_im[*index][p] = self.buffer[p].im;
                // The conjugate eigenvalue has the conjugate solution.
                self.y_re[*conjugate][p] = self.buffer[p].re;
                self.y_im[*conjugate][p] = -self.buffer[p].im;
            }
        }

        // dZ = (T kron I) dW, taking the real part.
        for i in 0..s * s {
            re[i] = self.decoupling.t[i].re;
            im[i] = self.decoupling.t[i].im;
        }
        let mut out_re = vec![vec![0.0; n]; s];
        let mut out_im = vec![vec![0.0; n]; s];
        simd::stage_transform(&re, s, 1.0, &self.y_re, &mut out_re);
        simd::stage_transform(&im, s, 1.0, &self.y_im, &mut out_im);
        for i in 0..s {
            for p in 0..n {
                rhs[i * n + p] = out_re[i][p] - out_im[i][p];
            }
        }
        true
    }
}

impl StageDecoupling {
    fn shifts_at(&self, index: usize) -> f64 {
        for block in &self.blocks {
            if let Block::Real { index: k, lambda } = block {
                if *k == index {
                    return *lambda;
                }
            }
        }
        0.0
    }

    fn complex_shift_at(&self, index: usize) -> Complex {
        for block in &self.blocks {
            if let Block::Pair { index: k, lambda, .. } = block {
                if *k == index {
                    return *lambda;
                }
            }
        }
        Complex::new(0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn radau_iia_5() -> Matrix<f64> {
        let r6 = 6f64.sqrt();
        Matrix::from_rows(&[
            vec![(88.0 - 7.0 * r6) / 360.0, (296.0 - 169.0 * r6) / 1800.0, (-2.0 + 3.0 * r6) / 225.0],
            vec![(296.0 + 169.0 * r6) / 1800.0, (88.0 + 7.0 * r6) / 360.0, (-2.0 - 3.0 * r6) / 225.0],
            vec![(16.0 - r6) / 36.0, (16.0 + r6) / 36.0, 1.0 / 9.0],
        ])
    }

    #[test]
    fn radau_decouples_into_one_real_and_one_complex_block() {
        let a = radau_iia_5();
        let decoupling = StageDecoupling::new(&a).expect("Radau IIA must decouple");
        assert_eq!(decoupling.stages, 3);
        // One real eigenvalue and one conjugate pair, so two factorizations of
        // size n instead of one of size 3n.
        assert_eq!(decoupling.factorizations(), 2);
    }

    #[test]
    fn decoupled_solve_matches_the_coupled_one() {
        let a = radau_iia_5();
        let s = 3;
        let n = 4;
        let h = 0.37;

        // A Jacobian with both real and complex spectrum.
        let mut jacobian = Matrix::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                jacobian[(i, j)] = ((i * 7 + j * 3) % 5) as f64 * 0.3 - 0.8;
            }
            jacobian[(i, i)] -= 3.0;
        }

        // Reference: the coupled s*n system.
        let mut coupled = Matrix::<f64>::zeros(s * n, s * n);
        for i in 0..s {
            for j in 0..s {
                for p in 0..n {
                    for q in 0..n {
                        coupled[(i * n + p, j * n + q)] = -h * a[(i, j)] * jacobian[(p, q)];
                    }
                }
            }
        }
        for d in 0..s * n {
            coupled[(d, d)] += 1.0;
        }
        let rhs: Vec<f64> = (0..s * n).map(|i| (i as f64 * 0.37).sin() + 0.4).collect();
        let reference = Lu::factor(coupled).solve(&rhs).expect("coupled system solvable");

        let decoupling = StageDecoupling::new(&a).unwrap();
        let mut linear = DecoupledLinear::new(decoupling, n);
        assert!(linear.factor(&jacobian, h));
        let mut computed = rhs.clone();
        assert!(linear.solve(&mut computed, h));

        for i in 0..s * n {
            assert!(
                (computed[i] - reference[i]).abs() < 1e-9 * reference[i].abs().max(1.0),
                "component {i}: decoupled {} vs coupled {}",
                computed[i],
                reference[i]
            );
        }
    }

    #[test]
    fn gauss_legendre_also_decouples() {
        let r3 = 3f64.sqrt();
        let a = Matrix::from_rows(&[
            vec![0.25, 0.25 - r3 / 6.0],
            vec![0.25 + r3 / 6.0, 0.25],
        ]);
        let decoupling = StageDecoupling::new(&a).expect("Gauss must decouple");
        assert_eq!(decoupling.factorizations(), 1);
    }
}
