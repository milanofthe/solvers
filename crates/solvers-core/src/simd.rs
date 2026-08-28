//! Vectorized kernels.
//!
//! Everything here is written in portable safe Rust in a shape that LLVM
//! reliably turns into packed instructions (SSE2/AVX on x86, NEON on ARM,
//! simd128 on wasm). No nightly features and no target specific intrinsics, so
//! the same code serves the native build and the browser build.
//!
//! The kernels that matter are the multi vector accumulations. A Runge-Kutta
//! stage update touches every stage vector once per component, and a fully
//! implicit stage system does the same across the whole `s * n` block. Blocking
//! over components and accumulating across stages inside the block keeps the
//! stage coefficients in registers and turns the update into a chain of fused
//! multiply-adds.

/// Components processed per block. Four doubles fill a 256 bit register and
/// split cleanly for narrower targets.
pub const LANES: usize = 4;

/// `y += a * x`
#[inline]
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    debug_assert_eq!(x.len(), y.len());
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi += a * *xi;
    }
}

/// `y = a * x`
#[inline]
pub fn scale_to(a: f64, x: &[f64], y: &mut [f64]) {
    debug_assert_eq!(x.len(), y.len());
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi = a * *xi;
    }
}

/// `x *= a`
#[inline]
pub fn scale(a: f64, x: &mut [f64]) {
    for xi in x.iter_mut() {
        *xi *= a;
    }
}

#[inline]
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), y.len());
    // Four partial sums so the additions are independent and can pipeline.
    let mut acc = [0.0f64; LANES];
    let mut chunks = x.chunks_exact(LANES).zip(y.chunks_exact(LANES));
    for (a, b) in &mut chunks {
        for lane in 0..LANES {
            acc[lane] += a[lane] * b[lane];
        }
    }
    let mut total = acc.iter().sum::<f64>();
    let rest = x.len() - x.len() % LANES;
    for i in rest..x.len() {
        total += x[i] * y[i];
    }
    total
}

/// Maximum absolute value.
#[inline]
pub fn norm_inf(x: &[f64]) -> f64 {
    x.iter().fold(0.0f64, |acc, v| acc.max(v.abs()))
}

/// Root mean square of `x` scaled component wise by `1 / scale`.
///
/// This is the norm the error controllers use, so it lives in one place.
#[inline]
pub fn weighted_rms(x: &[f64], scale: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), scale.len());
    if x.is_empty() {
        return 0.0;
    }
    let mut acc = [0.0f64; LANES];
    let mut chunks = x.chunks_exact(LANES).zip(scale.chunks_exact(LANES));
    for (a, s) in &mut chunks {
        for lane in 0..LANES {
            let r = a[lane] / s[lane];
            acc[lane] += r * r;
        }
    }
    let mut total = acc.iter().sum::<f64>();
    let rest = x.len() - x.len() % LANES;
    for i in rest..x.len() {
        let r = x[i] / scale[i];
        total += r * r;
    }
    (total / x.len() as f64).sqrt()
}

/// Fill `scale[i] = atol + rtol * max(|a[i]|, |b[i]|)`.
#[inline]
pub fn error_scale(atol: f64, rtol: f64, a: &[f64], b: &[f64], scale: &mut [f64]) {
    for i in 0..scale.len() {
        scale[i] = atol + rtol * a[i].abs().max(b[i].abs());
    }
}

/// `out = base + factor * sum_k coeffs[k] * vectors[k]`
///
/// The Runge-Kutta stage and solution update in one kernel. Blocking over
/// components lets the accumulator stay in registers across all stages.
#[inline]
pub fn combine(base: &[f64], factor: f64, coeffs: &[f64], vectors: &[Vec<f64>], out: &mut [f64]) {
    debug_assert_eq!(base.len(), out.len());
    debug_assert!(coeffs.len() <= vectors.len());
    let n = out.len();
    let blocks = n / LANES;

    for block in 0..blocks {
        let start = block * LANES;
        let mut acc = [0.0f64; LANES];
        for (k, &c) in coeffs.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            let v = &vectors[k][start..start + LANES];
            for lane in 0..LANES {
                acc[lane] += c * v[lane];
            }
        }
        for lane in 0..LANES {
            out[start + lane] = base[start + lane] + factor * acc[lane];
        }
    }

    for i in blocks * LANES..n {
        let mut acc = 0.0;
        for (k, &c) in coeffs.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            acc += c * vectors[k][i];
        }
        out[i] = base[i] + factor * acc;
    }
}

/// `out = factor * sum_k coeffs[k] * vectors[k]`
#[inline]
pub fn combine_into(factor: f64, coeffs: &[f64], vectors: &[Vec<f64>], out: &mut [f64]) {
    let n = out.len();
    let blocks = n / LANES;

    for block in 0..blocks {
        let start = block * LANES;
        let mut acc = [0.0f64; LANES];
        for (k, &c) in coeffs.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            let v = &vectors[k][start..start + LANES];
            for lane in 0..LANES {
                acc[lane] += c * v[lane];
            }
        }
        for lane in 0..LANES {
            out[start + lane] = factor * acc[lane];
        }
    }

    for i in blocks * LANES..n {
        let mut acc = 0.0;
        for (k, &c) in coeffs.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            acc += c * vectors[k][i];
        }
        out[i] = factor * acc;
    }
}

/// Apply a dense `s x s` coefficient matrix across `s` stage vectors of length
/// `n`: `out[i] = factor * sum_j m[i][j] * stages[j]`.
///
/// This is the transformation at the heart of a fully implicit stage system.
/// It is a small dense matrix times a tall block vector, so the profitable
/// blocking is over components with the whole stage dimension held inside.
#[inline]
pub fn stage_transform(m: &[f64], s: usize, factor: f64, stages: &[Vec<f64>], out: &mut [Vec<f64>]) {
    debug_assert_eq!(m.len(), s * s);
    debug_assert_eq!(stages.len(), s);
    debug_assert_eq!(out.len(), s);
    let n = stages[0].len();
    let blocks = n / LANES;

    for block in 0..blocks {
        let start = block * LANES;
        for i in 0..s {
            let mut acc = [0.0f64; LANES];
            for j in 0..s {
                let c = m[i * s + j];
                if c == 0.0 {
                    continue;
                }
                let v = &stages[j][start..start + LANES];
                for lane in 0..LANES {
                    acc[lane] += c * v[lane];
                }
            }
            let target = &mut out[i][start..start + LANES];
            for lane in 0..LANES {
                target[lane] = factor * acc[lane];
            }
        }
    }

    for idx in blocks * LANES..n {
        for i in 0..s {
            let mut acc = 0.0;
            for j in 0..s {
                acc += m[i * s + j] * stages[j][idx];
            }
            out[i][idx] = factor * acc;
        }
    }
}

/// Split a complex vector held as separate real and imaginary parts into the
/// interleaved layout the complex LU expects, and back.
#[inline]
pub fn interleave(re: &[f64], im: &[f64], out: &mut [crate::num::Complex]) {
    for i in 0..out.len() {
        out[i] = crate::num::Complex::new(re[i], im[i]);
    }
}

#[inline]
pub fn deinterleave(z: &[crate::num::Complex], re: &mut [f64], im: &mut [f64]) {
    for i in 0..z.len() {
        re[i] = z[i].re;
        im[i] = z[i].im;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combine_matches_the_naive_loop() {
        let n = 11;
        let base: Vec<f64> = (0..n).map(|i| i as f64 * 0.5).collect();
        let vectors: Vec<Vec<f64>> = (0..3)
            .map(|k| (0..n).map(|i| (i + k) as f64).collect())
            .collect();
        let coeffs = [0.25, 0.0, -1.5];
        let mut out = vec![0.0; n];
        combine(&base, 0.1, &coeffs, &vectors, &mut out);
        for i in 0..n {
            let expected = base[i]
                + 0.1 * (0.25 * vectors[0][i] + 0.0 * vectors[1][i] - 1.5 * vectors[2][i]);
            assert!((out[i] - expected).abs() < 1e-14);
        }
    }

    #[test]
    fn stage_transform_matches_the_naive_loop() {
        let (s, n) = (3, 7);
        let m: Vec<f64> = (0..s * s).map(|i| (i as f64) * 0.3 - 1.0).collect();
        let stages: Vec<Vec<f64>> = (0..s)
            .map(|k| (0..n).map(|i| (i as f64 + 1.0) * (k as f64 + 2.0)).collect())
            .collect();
        let mut out = vec![vec![0.0; n]; s];
        stage_transform(&m, s, 2.0, &stages, &mut out);
        for i in 0..s {
            for idx in 0..n {
                let mut expected = 0.0;
                for j in 0..s {
                    expected += m[i * s + j] * stages[j][idx];
                }
                assert!((out[i][idx] - 2.0 * expected).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn weighted_rms_is_the_usual_norm() {
        let x = [1.0, -2.0, 3.0, 4.0, 5.0];
        let s = [1.0; 5];
        let expected = (x.iter().map(|v| v * v).sum::<f64>() / 5.0).sqrt();
        assert!((weighted_rms(&x, &s) - expected).abs() < 1e-14);
    }
}
