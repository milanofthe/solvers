//! Exact rational arithmetic on i128, with checked operations.
//!
//! Butcher tableaux are published as exact fractions. Keeping them exact lets
//! the order conditions be verified as identities instead of being measured
//! against a tolerance.

use std::cmp::Ordering;

#[derive(Copy, Clone, Debug)]
pub struct Rational {
    num: i128,
    den: i128,
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    if a == 0 {
        1
    } else {
        a
    }
}

impl Rational {
    pub const ZERO: Rational = Rational { num: 0, den: 1 };
    pub const ONE: Rational = Rational { num: 1, den: 1 };

    /// Build a normalized rational, `None` if the denominator is zero.
    pub fn new(num: i128, den: i128) -> Option<Rational> {
        if den == 0 {
            return None;
        }
        let sign = if den < 0 { -1 } else { 1 };
        let g = gcd(num, den);
        Some(Rational {
            num: sign * (num / g),
            den: sign * (den / g),
        })
    }

    pub fn integer(n: i128) -> Rational {
        Rational { num: n, den: 1 }
    }

    pub fn numerator(self) -> i128 {
        self.num
    }

    pub fn denominator(self) -> i128 {
        self.den
    }

    pub fn is_zero(self) -> bool {
        self.num == 0
    }

    pub fn to_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }

    pub fn checked_add(self, rhs: Rational) -> Option<Rational> {
        let g = gcd(self.den, rhs.den);
        let lcm = self.den.checked_div(g)?.checked_mul(rhs.den)?;
        let a = self.num.checked_mul(lcm / self.den)?;
        let b = rhs.num.checked_mul(lcm / rhs.den)?;
        Rational::new(a.checked_add(b)?, lcm)
    }

    pub fn checked_sub(self, rhs: Rational) -> Option<Rational> {
        self.checked_add(rhs.neg())
    }

    pub fn checked_mul(self, rhs: Rational) -> Option<Rational> {
        // Cross reduce first to keep the intermediates small.
        let g1 = gcd(self.num, rhs.den);
        let g2 = gcd(rhs.num, self.den);
        let num = (self.num / g1).checked_mul(rhs.num / g2)?;
        let den = (self.den / g2).checked_mul(rhs.den / g1)?;
        Rational::new(num, den)
    }

    pub fn checked_div(self, rhs: Rational) -> Option<Rational> {
        if rhs.num == 0 {
            return None;
        }
        self.checked_mul(Rational {
            num: rhs.den,
            den: rhs.num,
        })
    }

    pub fn neg(self) -> Rational {
        Rational {
            num: -self.num,
            den: self.den,
        }
    }

    pub fn checked_powi(self, n: i32) -> Option<Rational> {
        let mut acc = Rational::ONE;
        let base = if n < 0 {
            Rational::ONE.checked_div(self)?
        } else {
            self
        };
        for _ in 0..n.unsigned_abs() {
            acc = acc.checked_mul(base)?;
        }
        Some(acc)
    }

    /// Rationalize a float via continued fractions.
    ///
    /// Only accepted when the result reproduces the input bit for bit, so a
    /// literal `0.25` in JSON stays exact while a truncated decimal does not
    /// silently turn into a wrong fraction.
    pub fn from_f64_exact(v: f64) -> Option<Rational> {
        if !v.is_finite() {
            return None;
        }
        if v == 0.0 {
            return Some(Rational::ZERO);
        }
        const MAX_DEN: i128 = 1_000_000_000_000;
        let (mut h0, mut h1) = (0i128, 1i128);
        let (mut k0, mut k1) = (1i128, 0i128);
        let mut x = v;
        for _ in 0..64 {
            let a = x.floor();
            if a.abs() > 1e18 {
                return None;
            }
            let ai = a as i128;
            let h2 = ai.checked_mul(h1)?.checked_add(h0)?;
            let k2 = ai.checked_mul(k1)?.checked_add(k0)?;
            if k2.abs() > MAX_DEN {
                break;
            }
            h0 = h1;
            h1 = h2;
            k0 = k1;
            k1 = k2;
            let candidate = Rational::new(h1, k1)?;
            if candidate.to_f64() == v {
                return Some(candidate);
            }
            let frac = x - a;
            if frac == 0.0 {
                break;
            }
            x = 1.0 / frac;
        }
        None
    }
}

impl PartialEq for Rational {
    fn eq(&self, other: &Self) -> bool {
        self.num == other.num && self.den == other.den
    }
}

impl Eq for Rational {}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        // Denominators are positive by construction.
        match (self.num as i128).checked_mul(other.den) {
            Some(lhs) => match (other.num as i128).checked_mul(self.den) {
                Some(rhs) => lhs.cmp(&rhs),
                None => self.to_f64().total_cmp(&other.to_f64()),
            },
            None => self.to_f64().total_cmp(&other.to_f64()),
        }
    }
}

impl std::fmt::Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.den == 1 {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}
