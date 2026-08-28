//! Deserialization helpers for coefficients in method files.
//!
//! A coefficient may be written as a JSON number (`0.25`), as an arithmetic
//! expression (`"(2 - sqrt(2))/4"`), or, for multistep families, as the marker
//! `"free"` meaning "determine me from the order conditions".

use crate::num::{expr, Coeff};
use serde::de::{self, Deserializer, Visitor};
use serde::{Deserialize, Serialize, Serializer};
use std::fmt;

/// A coefficient with a value.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CoeffValue(pub Coeff);

impl CoeffValue {
    pub fn get(self) -> Coeff {
        self.0
    }
    pub fn value(self) -> f64 {
        self.0.value()
    }
}

impl From<Coeff> for CoeffValue {
    fn from(c: Coeff) -> Self {
        CoeffValue(c)
    }
}

struct CoeffVisitor;

impl<'de> Visitor<'de> for CoeffVisitor {
    type Value = Coeff;

    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("a number or an arithmetic expression string")
    }

    fn visit_i64<E: de::Error>(self, v: i64) -> Result<Coeff, E> {
        Ok(Coeff::rational(v as i128, 1))
    }

    fn visit_u64<E: de::Error>(self, v: u64) -> Result<Coeff, E> {
        Ok(Coeff::rational(v as i128, 1))
    }

    fn visit_f64<E: de::Error>(self, v: f64) -> Result<Coeff, E> {
        Ok(Coeff::from_f64_rationalized(v))
    }

    fn visit_str<E: de::Error>(self, v: &str) -> Result<Coeff, E> {
        expr::eval(v).map_err(|e| E::custom(format!("invalid coefficient {v:?}: {e}")))
    }
}

impl<'de> Deserialize<'de> for CoeffValue {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        d.deserialize_any(CoeffVisitor).map(CoeffValue)
    }
}

impl Serialize for CoeffValue {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self.0 {
            Coeff::Exact(r) => s.serialize_str(&r.to_string()),
            Coeff::Approx(v) => s.serialize_f64(v),
        }
    }
}

/// Either a fixed coefficient or an unknown to be solved for.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Slot {
    Free,
    Fixed(Coeff),
}

impl Slot {
    pub fn is_free(self) -> bool {
        matches!(self, Slot::Free)
    }
    pub fn fixed_value(self) -> f64 {
        match self {
            Slot::Free => 0.0,
            Slot::Fixed(c) => c.value(),
        }
    }
}

struct SlotVisitor;

impl<'de> Visitor<'de> for SlotVisitor {
    type Value = Slot;

    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("a number, an expression string, or \"free\"")
    }

    fn visit_i64<E: de::Error>(self, v: i64) -> Result<Slot, E> {
        Ok(Slot::Fixed(Coeff::rational(v as i128, 1)))
    }

    fn visit_u64<E: de::Error>(self, v: u64) -> Result<Slot, E> {
        Ok(Slot::Fixed(Coeff::rational(v as i128, 1)))
    }

    fn visit_f64<E: de::Error>(self, v: f64) -> Result<Slot, E> {
        Ok(Slot::Fixed(Coeff::from_f64_rationalized(v)))
    }

    fn visit_unit<E: de::Error>(self) -> Result<Slot, E> {
        Ok(Slot::Free)
    }

    fn visit_none<E: de::Error>(self) -> Result<Slot, E> {
        Ok(Slot::Free)
    }

    fn visit_str<E: de::Error>(self, v: &str) -> Result<Slot, E> {
        let trimmed = v.trim();
        if trimmed.eq_ignore_ascii_case("free") || trimmed == "?" {
            return Ok(Slot::Free);
        }
        expr::eval(trimmed)
            .map(Slot::Fixed)
            .map_err(|e| E::custom(format!("invalid coefficient {v:?}: {e}")))
    }
}

impl<'de> Deserialize<'de> for Slot {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        d.deserialize_any(SlotVisitor)
    }
}

impl Serialize for Slot {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            Slot::Free => s.serialize_str("free"),
            Slot::Fixed(c) => CoeffValue(*c).serialize(s),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_numbers_and_expressions() {
        let v: Vec<CoeffValue> = serde_json::from_str(r#"[0.25, "1/3", 2, "(2-sqrt(2))/4"]"#).unwrap();
        assert!(v[0].get().is_exact());
        assert!(v[1].get().is_exact());
        assert_eq!(v[2].value(), 2.0);
        assert!(!v[3].get().is_exact());
    }

    #[test]
    fn reads_free_slots() {
        let v: Vec<Slot> = serde_json::from_str(r#"["free", 0, null, "1/2"]"#).unwrap();
        assert!(v[0].is_free());
        assert!(!v[1].is_free());
        assert!(v[2].is_free());
        assert_eq!(v[3].fixed_value(), 0.5);
    }
}
