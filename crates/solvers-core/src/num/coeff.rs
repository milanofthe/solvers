//! Tableau coefficients that stay exact for as long as they can.
//!
//! `Coeff` is a rational that degrades to a float on overflow or as soon as an
//! irrational operation such as `sqrt` is involved. Order conditions are then
//! checked exactly for the (many) purely rational methods and numerically for
//! the rest, without the caller having to care which case it is in.

use super::{Field, Rational};
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Copy, Clone, Debug)]
pub enum Coeff {
    Exact(Rational),
    Approx(f64),
}

impl Coeff {
    pub fn rational(num: i128, den: i128) -> Coeff {
        match Rational::new(num, den) {
            Some(r) => Coeff::Exact(r),
            None => Coeff::Approx(f64::NAN),
        }
    }

    pub fn value(self) -> f64 {
        match self {
            Coeff::Exact(r) => r.to_f64(),
            Coeff::Approx(v) => v,
        }
    }

    pub fn is_exact(self) -> bool {
        matches!(self, Coeff::Exact(_))
    }

    pub fn exact(self) -> Option<Rational> {
        match self {
            Coeff::Exact(r) => Some(r),
            Coeff::Approx(_) => None,
        }
    }

    /// Adopt a float, keeping it exact when it rationalizes without loss.
    pub fn from_f64_rationalized(v: f64) -> Coeff {
        match Rational::from_f64_exact(v) {
            Some(r) => Coeff::Exact(r),
            None => Coeff::Approx(v),
        }
    }

    fn binary(
        self,
        rhs: Coeff,
        exact: impl Fn(Rational, Rational) -> Option<Rational>,
        approx: impl Fn(f64, f64) -> f64,
    ) -> Coeff {
        if let (Coeff::Exact(a), Coeff::Exact(b)) = (self, rhs) {
            if let Some(r) = exact(a, b) {
                return Coeff::Exact(r);
            }
        }
        Coeff::Approx(approx(self.value(), rhs.value()))
    }

    pub fn powi(self, n: i32) -> Coeff {
        if let Coeff::Exact(r) = self {
            if let Some(p) = r.checked_powi(n) {
                return Coeff::Exact(p);
            }
        }
        Coeff::Approx(self.value().powi(n))
    }

    pub fn sqrt(self) -> Coeff {
        // Perfect squares of rationals stay exact.
        if let Coeff::Exact(r) = self {
            let (n, d) = (r.numerator(), r.denominator());
            if n >= 0 {
                let sn = (n as f64).sqrt().round() as i128;
                let sd = (d as f64).sqrt().round() as i128;
                if sn * sn == n && sd * sd == d {
                    if let Some(q) = Rational::new(sn, sd) {
                        return Coeff::Exact(q);
                    }
                }
            }
        }
        Coeff::Approx(self.value().sqrt())
    }
}

impl Add for Coeff {
    type Output = Coeff;
    fn add(self, o: Coeff) -> Coeff {
        self.binary(o, |a, b| a.checked_add(b), |a, b| a + b)
    }
}

impl Sub for Coeff {
    type Output = Coeff;
    fn sub(self, o: Coeff) -> Coeff {
        self.binary(o, |a, b| a.checked_sub(b), |a, b| a - b)
    }
}

impl Mul for Coeff {
    type Output = Coeff;
    fn mul(self, o: Coeff) -> Coeff {
        self.binary(o, |a, b| a.checked_mul(b), |a, b| a * b)
    }
}

impl Div for Coeff {
    type Output = Coeff;
    fn div(self, o: Coeff) -> Coeff {
        self.binary(o, |a, b| a.checked_div(b), |a, b| a / b)
    }
}

impl Neg for Coeff {
    type Output = Coeff;
    fn neg(self) -> Coeff {
        match self {
            Coeff::Exact(r) => Coeff::Exact(r.neg()),
            Coeff::Approx(v) => Coeff::Approx(-v),
        }
    }
}

impl PartialEq for Coeff {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Coeff::Exact(a), Coeff::Exact(b)) => a == b,
            _ => self.value() == other.value(),
        }
    }
}

impl Field for Coeff {
    fn zero() -> Self {
        Coeff::Exact(Rational::ZERO)
    }
    fn one() -> Self {
        Coeff::Exact(Rational::ONE)
    }
    fn from_f64(v: f64) -> Self {
        Coeff::Approx(v)
    }
    fn from_i64(v: i64) -> Self {
        Coeff::Exact(Rational::integer(v as i128))
    }
    fn magnitude(self) -> f64 {
        self.value().abs()
    }
    fn is_zero(self) -> bool {
        match self {
            Coeff::Exact(r) => r.is_zero(),
            Coeff::Approx(v) => v == 0.0,
        }
    }
}

impl std::fmt::Display for Coeff {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Coeff::Exact(r) => write!(f, "{}", r),
            Coeff::Approx(v) => write!(f, "{}", v),
        }
    }
}
