//! A very small arithmetic expression evaluator for method parameter files.
//!
//! Published tableaux mix plain fractions with closed form irrationals, so the
//! JSON files accept strings such as `"2012122486997/3467029789466"`,
//! `"(2 - sqrt(2))/4"` or `"1/(2 + sqrt(2))"`. Everything that can stay a
//! rational does, the rest falls back to floating point.

use super::{Coeff, Rational};

#[derive(Debug, Clone, PartialEq)]
pub struct ExprError {
    pub message: String,
    pub position: usize,
}

impl std::fmt::Display for ExprError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at position {}", self.message, self.position)
    }
}

impl std::error::Error for ExprError {}

/// Evaluate an arithmetic expression to a coefficient.
pub fn eval(input: &str) -> Result<Coeff, ExprError> {
    let chars: Vec<char> = input.chars().collect();
    let mut p = Parser { chars, pos: 0 };
    p.skip_ws();
    let value = p.expr()?;
    p.skip_ws();
    if p.pos < p.chars.len() {
        return Err(p.error("unexpected trailing input"));
    }
    Ok(value)
}

struct Parser {
    chars: Vec<char>,
    pos: usize,
}

impl Parser {
    fn error(&self, message: &str) -> ExprError {
        ExprError {
            message: message.to_string(),
            position: self.pos,
        }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(c) if c.is_whitespace()) {
            self.pos += 1;
        }
    }

    fn eat(&mut self, c: char) -> bool {
        self.skip_ws();
        if self.peek() == Some(c) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn expr(&mut self) -> Result<Coeff, ExprError> {
        let mut lhs = self.term()?;
        loop {
            self.skip_ws();
            match self.peek() {
                Some('+') => {
                    self.pos += 1;
                    lhs = lhs + self.term()?;
                }
                Some('-') => {
                    self.pos += 1;
                    lhs = lhs - self.term()?;
                }
                _ => return Ok(lhs),
            }
        }
    }

    fn term(&mut self) -> Result<Coeff, ExprError> {
        let mut lhs = self.unary()?;
        loop {
            self.skip_ws();
            match self.peek() {
                Some('*') => {
                    self.pos += 1;
                    lhs = lhs * self.unary()?;
                }
                Some('/') => {
                    self.pos += 1;
                    lhs = lhs / self.unary()?;
                }
                _ => return Ok(lhs),
            }
        }
    }

    fn unary(&mut self) -> Result<Coeff, ExprError> {
        self.skip_ws();
        match self.peek() {
            Some('-') => {
                self.pos += 1;
                Ok(-self.unary()?)
            }
            Some('+') => {
                self.pos += 1;
                self.unary()
            }
            _ => self.power(),
        }
    }

    fn power(&mut self) -> Result<Coeff, ExprError> {
        let base = self.atom()?;
        self.skip_ws();
        if self.peek() == Some('^') {
            self.pos += 1;
            let exponent = self.unary()?;
            let e = exponent.value();
            if e.fract() == 0.0 && e.abs() <= 64.0 {
                return Ok(base.powi(e as i32));
            }
            return Ok(Coeff::Approx(base.value().powf(e)));
        }
        Ok(base)
    }

    fn atom(&mut self) -> Result<Coeff, ExprError> {
        self.skip_ws();
        match self.peek() {
            None => Err(self.error("unexpected end of expression")),
            Some('(') => {
                self.pos += 1;
                let inner = self.expr()?;
                if !self.eat(')') {
                    return Err(self.error("missing closing parenthesis"));
                }
                Ok(inner)
            }
            Some(c) if c.is_ascii_digit() || c == '.' => self.number(),
            Some(c) if c.is_alphabetic() || c == '_' => self.identifier(),
            Some(_) => Err(self.error("unexpected character")),
        }
    }

    fn identifier(&mut self) -> Result<Coeff, ExprError> {
        let start = self.pos;
        while matches!(self.peek(), Some(c) if c.is_alphanumeric() || c == '_') {
            self.pos += 1;
        }
        let name: String = self.chars[start..self.pos].iter().collect();
        let name = name.to_ascii_lowercase();

        // Constants take no argument list.
        match name.as_str() {
            "pi" => return Ok(Coeff::Approx(std::f64::consts::PI)),
            "e" => return Ok(Coeff::Approx(std::f64::consts::E)),
            _ => {}
        }

        if !self.eat('(') {
            return Err(self.error("unknown constant"));
        }
        let arg = self.expr()?;
        if !self.eat(')') {
            return Err(self.error("missing closing parenthesis"));
        }

        let v = arg.value();
        Ok(match name.as_str() {
            "sqrt" => arg.sqrt(),
            "cbrt" => Coeff::Approx(v.cbrt()),
            "exp" => Coeff::Approx(v.exp()),
            "ln" | "log" => Coeff::Approx(v.ln()),
            "log10" => Coeff::Approx(v.log10()),
            "sin" => Coeff::Approx(v.sin()),
            "cos" => Coeff::Approx(v.cos()),
            "tan" => Coeff::Approx(v.tan()),
            "atan" => Coeff::Approx(v.atan()),
            "abs" => {
                if v < 0.0 {
                    -arg
                } else {
                    arg
                }
            }
            _ => return Err(self.error("unknown function")),
        })
    }

    /// Decimal literals are read exactly: `0.25` becomes `1/4`, not a float.
    fn number(&mut self) -> Result<Coeff, ExprError> {
        let start = self.pos;
        let mut digits = String::new();
        let mut scale: i32 = 0;
        let mut seen_dot = false;

        while let Some(c) = self.peek() {
            if c.is_ascii_digit() {
                digits.push(c);
                if seen_dot {
                    scale += 1;
                }
                self.pos += 1;
            } else if c == '.' && !seen_dot {
                seen_dot = true;
                self.pos += 1;
            } else if c == '_' {
                self.pos += 1;
            } else {
                break;
            }
        }
        if digits.is_empty() {
            self.pos = start;
            return Err(self.error("malformed number"));
        }

        let mut exponent: i32 = 0;
        if matches!(self.peek(), Some('e') | Some('E')) {
            // Only consume as an exponent when it really looks like one.
            let save = self.pos;
            self.pos += 1;
            let mut sign = 1i32;
            match self.peek() {
                Some('+') => self.pos += 1,
                Some('-') => {
                    sign = -1;
                    self.pos += 1;
                }
                _ => {}
            }
            let mut exp_digits = String::new();
            while matches!(self.peek(), Some(c) if c.is_ascii_digit()) {
                exp_digits.push(self.chars[self.pos]);
                self.pos += 1;
            }
            if exp_digits.is_empty() {
                self.pos = save;
            } else {
                exponent = sign * exp_digits.parse::<i32>().unwrap_or(0);
            }
        }

        let text: String = self.chars[start..self.pos].iter().collect();
        let net = exponent - scale;

        if let Ok(mantissa) = digits.parse::<i128>() {
            let base = Rational::integer(mantissa);
            let ten = Rational::integer(10);
            let scaled = if net >= 0 {
                ten.checked_powi(net).and_then(|p| base.checked_mul(p))
            } else {
                ten.checked_powi(-net).and_then(|p| base.checked_div(p))
            };
            if let Some(r) = scaled {
                return Ok(Coeff::Exact(r));
            }
        }

        text.parse::<f64>()
            .map(Coeff::Approx)
            .map_err(|_| self.error("malformed number"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(s: &str) -> f64 {
        eval(s).unwrap().value()
    }

    #[test]
    fn fractions_stay_exact() {
        let c = eval("1/4").unwrap();
        assert!(c.is_exact());
        assert_eq!(c.value(), 0.25);
    }

    #[test]
    fn decimals_are_rationalized() {
        let c = eval("0.25").unwrap();
        assert!(c.is_exact());
        assert_eq!(c.exact().unwrap().denominator(), 4);
    }

    #[test]
    fn long_fractions_survive() {
        let c = eval("2012122486997/3467029789466").unwrap();
        assert!(c.is_exact());
    }

    #[test]
    fn irrationals_degrade() {
        let c = eval("(2 - sqrt(2))/4").unwrap();
        assert!(!c.is_exact());
        assert!((c.value() - (2.0 - 2f64.sqrt()) / 4.0).abs() < 1e-15);
    }

    #[test]
    fn perfect_squares_stay_exact() {
        assert!(eval("sqrt(9/4)").unwrap().is_exact());
    }

    #[test]
    fn precedence_and_signs() {
        assert_eq!(approx("1 - 2*3"), -5.0);
        assert_eq!(approx("2^3^2"), 512.0);
        assert_eq!(approx("-3/28"), -3.0 / 28.0);
        assert_eq!(approx("1.5e2"), 150.0);
    }
}
