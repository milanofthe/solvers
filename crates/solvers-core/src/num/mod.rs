//! Scalar types.
//!
//! The framework is generic over a small `Field` abstraction so that the same
//! linear algebra and the same tableau evaluation code can run on real numbers
//! (time stepping), complex numbers (stability functions) and exact rationals
//! (order condition verification).

mod coeff;
mod complex;
pub mod expr;
mod rational;

pub use coeff::Coeff;
pub use complex::Complex;
pub use rational::Rational;

use std::ops::{Add, Div, Mul, Neg, Sub};

/// Minimal field abstraction: everything the dense linear algebra needs.
pub trait Field:
    Copy
    + Clone
    + PartialEq
    + std::fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_f64(v: f64) -> Self;

    /// Non negative size used for pivot selection and zero tests.
    fn magnitude(self) -> f64;

    fn is_zero(self) -> bool {
        self.magnitude() == 0.0
    }

    fn from_i64(v: i64) -> Self {
        Self::from_f64(v as f64)
    }

    fn powi(self, n: u32) -> Self {
        let mut acc = Self::one();
        for _ in 0..n {
            acc = acc * self;
        }
        acc
    }
}

impl Field for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_f64(v: f64) -> Self {
        v
    }
    fn magnitude(self) -> f64 {
        self.abs()
    }
}
