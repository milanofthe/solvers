//! Complex numbers.
//!
//! Only what the stability analysis needs: a field implementation so the dense
//! LU solver can be reused for `R(z) = 1 + z b^T (I - zA)^{-1} e`.

use super::Field;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub const I: Complex = Complex { re: 0.0, im: 1.0 };

    pub fn new(re: f64, im: f64) -> Complex {
        Complex { re, im }
    }

    pub fn real(re: f64) -> Complex {
        Complex { re, im: 0.0 }
    }

    /// `exp(i * theta)` on the unit circle.
    pub fn expi(theta: f64) -> Complex {
        Complex::new(theta.cos(), theta.sin())
    }

    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    pub fn conj(self) -> Complex {
        Complex::new(self.re, -self.im)
    }

    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }

    pub fn exp(self) -> Complex {
        let m = self.re.exp();
        Complex::new(m * self.im.cos(), m * self.im.sin())
    }

    pub fn sqrt(self) -> Complex {
        let r = self.abs();
        if r == 0.0 {
            return Complex::real(0.0);
        }
        let re = ((r + self.re) * 0.5).max(0.0).sqrt();
        let im = ((r - self.re) * 0.5).max(0.0).sqrt();
        Complex::new(re, if self.im < 0.0 { -im } else { im })
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, o: Complex) -> Complex {
        Complex::new(self.re + o.re, self.im + o.im)
    }
}

impl Sub for Complex {
    type Output = Complex;
    fn sub(self, o: Complex) -> Complex {
        Complex::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for Complex {
    type Output = Complex;
    fn mul(self, o: Complex) -> Complex {
        Complex::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

impl Div for Complex {
    type Output = Complex;
    fn div(self, o: Complex) -> Complex {
        // Smith's algorithm, avoids overflow when one part dominates.
        if o.re.abs() >= o.im.abs() {
            let r = o.im / o.re;
            let d = o.re + o.im * r;
            Complex::new((self.re + self.im * r) / d, (self.im - self.re * r) / d)
        } else {
            let r = o.re / o.im;
            let d = o.re * r + o.im;
            Complex::new((self.re * r + self.im) / d, (self.im * r - self.re) / d)
        }
    }
}

impl Neg for Complex {
    type Output = Complex;
    fn neg(self) -> Complex {
        Complex::new(-self.re, -self.im)
    }
}

impl Field for Complex {
    fn zero() -> Self {
        Complex::new(0.0, 0.0)
    }
    fn one() -> Self {
        Complex::new(1.0, 0.0)
    }
    fn from_f64(v: f64) -> Self {
        Complex::real(v)
    }
    fn magnitude(self) -> f64 {
        self.abs()
    }
}

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im >= 0.0 {
            write!(f, "{}+{}i", self.re, self.im)
        } else {
            write!(f, "{}{}i", self.re, self.im)
        }
    }
}
