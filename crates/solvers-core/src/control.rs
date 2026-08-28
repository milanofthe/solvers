//! Step size control.
//!
//! Every adaptive method in this crate uses the same controller, written in
//! Soederlind's general digital filter form
//!
//! ```text
//! h_{n+1} / h_n = S * (1/e_n)^(b1/k) * (1/e_{n-1})^(b2/k) * (1/e_{n-2})^(b3/k)
//!                   * (h_n / h_{n-1})^(-a2) * (h_{n-1} / h_{n-2})^(-a3)
//! ```
//!
//! where `e` is the error estimate already normalized so that `e <= 1` means
//! accept, and `k` is the order the error estimate converges with. The named
//! controllers in the literature are just parameter sets of this one formula,
//! which is why they are data here and not separate implementations.
//!
//! References
//! ----------
//! * G. Soederlind, "Digital filters in adaptive time-stepping",
//!   ACM TOMS 29(1), 2003, doi:10.1145/641876.641877
//! * K. Gustafsson, "Control theoretic techniques for stepsize selection in
//!   implicit Runge-Kutta methods", ACM TOMS 20(4), 1994,
//!   doi:10.1145/198429.198437
//! * E. Hairer, S. P. Noersett, G. Wanner, "Solving Ordinary Differential
//!   Equations I", 2nd ed., Springer 1993, doi:10.1007/978-3-540-78862-1

use serde::{Deserialize, Serialize};

/// Named parameter sets. `Custom` keeps whatever the caller filled in.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ControllerPreset {
    /// Deadbeat integral control, the classical `h * (1/e)^(1/k)`.
    I,
    /// Soederlind PI.3.3.3, a mild proportional integral filter.
    Pi3333,
    /// Soederlind PI.4.2, the common choice for explicit Runge-Kutta pairs.
    Pi4020,
    /// Hairer and Wanner's PI controller as used in DOPRI5.
    PiHairer,
    /// Soederlind H211PI, smoother than PI at the cost of response speed.
    H211Pi,
    /// Soederlind H211b with b = 4, includes step size ratio damping.
    H211b,
    /// Soederlind H312PID.
    H312Pid,
    /// Soederlind H312b with b = 8.
    H312b,
    /// Classical PID parameters.
    Pid,
    /// Gustafsson's predictive controller, the default of Radau5.
    Gustafsson,
    Custom,
}

impl ControllerPreset {
    pub fn from_name(name: &str) -> Option<ControllerPreset> {
        Some(match name.trim().to_ascii_lowercase().replace(['-', '_', '.'], "").as_str() {
            "i" | "integral" | "elementary" | "standard" => ControllerPreset::I,
            "pi3333" => ControllerPreset::Pi3333,
            "pi4020" | "pi42" => ControllerPreset::Pi4020,
            "pihairer" | "hairer" => ControllerPreset::PiHairer,
            "h211pi" => ControllerPreset::H211Pi,
            "h211b" => ControllerPreset::H211b,
            "h312pid" => ControllerPreset::H312Pid,
            "h312b" => ControllerPreset::H312b,
            "pid" => ControllerPreset::Pid,
            "gustafsson" | "predictive" => ControllerPreset::Gustafsson,
            "custom" => ControllerPreset::Custom,
            _ => return None,
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            ControllerPreset::I => "i",
            ControllerPreset::Pi3333 => "pi3333",
            ControllerPreset::Pi4020 => "pi4020",
            ControllerPreset::PiHairer => "pi_hairer",
            ControllerPreset::H211Pi => "h211pi",
            ControllerPreset::H211b => "h211b",
            ControllerPreset::H312Pid => "h312pid",
            ControllerPreset::H312b => "h312b",
            ControllerPreset::Pid => "pid",
            ControllerPreset::Gustafsson => "gustafsson",
            ControllerPreset::Custom => "custom",
        }
    }

    pub fn all() -> &'static [ControllerPreset] {
        use ControllerPreset::*;
        &[I, Pi3333, Pi4020, PiHairer, H211Pi, H211b, H312Pid, H312b, Pid, Gustafsson]
    }

    /// Error exponents `b` and step ratio exponents `a`.
    fn gains(self) -> ([f64; 3], [f64; 2]) {
        match self {
            ControllerPreset::I | ControllerPreset::Custom => ([1.0, 0.0, 0.0], [0.0, 0.0]),
            ControllerPreset::Pi3333 => ([2.0 / 3.0, -1.0 / 3.0, 0.0], [0.0, 0.0]),
            ControllerPreset::Pi4020 => ([3.0 / 5.0, -1.0 / 5.0, 0.0], [0.0, 0.0]),
            ControllerPreset::PiHairer => ([0.7, -0.4, 0.0], [0.0, 0.0]),
            ControllerPreset::H211Pi => ([1.0 / 6.0, 1.0 / 6.0, 0.0], [0.0, 0.0]),
            ControllerPreset::H211b => ([1.0 / 4.0, 1.0 / 4.0, 0.0], [1.0 / 4.0, 0.0]),
            ControllerPreset::H312Pid => ([1.0 / 18.0, 1.0 / 9.0, 1.0 / 18.0], [0.0, 0.0]),
            ControllerPreset::H312b => ([1.0 / 8.0, 2.0 / 8.0, 1.0 / 8.0], [3.0 / 8.0, 1.0 / 8.0]),
            ControllerPreset::Pid => ([0.58, -0.21, 0.10], [0.0, 0.0]),
            // h_{n+1} = h_n * (h_n / h_{n-1}) * (e_{n-1} / e_n^2)^(1/k)
            ControllerPreset::Gustafsson => ([2.0, -1.0, 0.0], [-1.0, 0.0]),
        }
    }
}

/// Controller parameters.
#[derive(Copy, Clone, Debug, Deserialize, Serialize)]
pub struct ControllerConfig {
    pub preset: ControllerPreset,
    /// Error exponents, only used when `preset` is `Custom`.
    pub beta: [f64; 3],
    /// Step size ratio exponents, only used when `preset` is `Custom`.
    pub alpha: [f64; 2],
    /// Safety factor applied to every proposal.
    pub safety: f64,
    pub min_scale: f64,
    pub max_scale: f64,
    /// Cap on the growth in the step directly after a rejection.
    pub max_scale_after_reject: f64,
    /// Soederlind's smooth limiter width. `None` uses hard clamping only.
    pub limiter: Option<f64>,
    /// Fall back to plain integral control while steps are being rejected.
    pub integral_on_reject: bool,
}

impl Default for ControllerConfig {
    fn default() -> Self {
        ControllerConfig {
            preset: ControllerPreset::Pi4020,
            beta: [1.0, 0.0, 0.0],
            alpha: [0.0, 0.0],
            safety: 0.9,
            min_scale: 0.2,
            max_scale: 10.0,
            max_scale_after_reject: 1.0,
            limiter: None,
            integral_on_reject: true,
        }
    }
}

impl ControllerConfig {
    pub fn preset(preset: ControllerPreset) -> ControllerConfig {
        let (beta, alpha) = preset.gains();
        ControllerConfig {
            preset,
            beta,
            alpha,
            ..Default::default()
        }
    }

    pub fn from_name(name: &str) -> Option<ControllerConfig> {
        ControllerPreset::from_name(name).map(ControllerConfig::preset)
    }

    fn gains(&self) -> ([f64; 3], [f64; 2]) {
        if self.preset == ControllerPreset::Custom {
            (self.beta, self.alpha)
        } else {
            self.preset.gains()
        }
    }
}

/// What the controller decided about a step.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Decision {
    pub accepted: bool,
    /// Multiplier for the next step size.
    pub scale: f64,
}

/// A controller with its filter memory.
#[derive(Clone, Debug)]
pub struct Controller {
    pub config: ControllerConfig,
    /// Previous error norms, most recent first.
    errors: [f64; 2],
    /// Previous step size ratios, most recent first.
    ratios: [f64; 2],
    history: usize,
    rejected_in_a_row: u32,
}

impl Controller {
    pub fn new(config: ControllerConfig) -> Controller {
        Controller {
            config,
            errors: [1.0, 1.0],
            ratios: [1.0, 1.0],
            history: 0,
            rejected_in_a_row: 0,
        }
    }

    pub fn reset(&mut self) {
        self.errors = [1.0, 1.0];
        self.ratios = [1.0, 1.0];
        self.history = 0;
        self.rejected_in_a_row = 0;
    }

    pub fn rejections(&self) -> u32 {
        self.rejected_in_a_row
    }

    /// Judge a step from its normalized error and return the step size factor.
    ///
    /// `order` is the order the error estimate converges with, so the exponent
    /// denominator is `order + 1`.
    pub fn control(&mut self, error: f64, step_ratio: f64, order: usize) -> Decision {
        let k = (order.max(1) + 1) as f64;
        let accepted = error <= 1.0;

        // A vanishing error means the step is unconstrained by accuracy.
        if !error.is_finite() {
            self.rejected_in_a_row += 1;
            return Decision {
                accepted: false,
                scale: self.config.min_scale,
            };
        }
        let e0 = error.max(1e-10);

        let use_integral =
            self.history == 0 || (!accepted && self.config.integral_on_reject) || self.rejected_in_a_row > 0;

        let (beta, alpha) = self.config.gains();
        let mut raw = if use_integral {
            e0.powf(-1.0 / k)
        } else {
            let e1 = self.errors[0].max(1e-10);
            let e2 = self.errors[1].max(1e-10);
            let mut f = e0.powf(-beta[0] / k) * e1.powf(-beta[1] / k);
            if self.history >= 2 {
                f *= e2.powf(-beta[2] / k);
                f *= self.ratios[1].powf(-alpha[1]);
            }
            f *= self.ratios[0].powf(-alpha[0]);
            f
        };

        raw *= self.config.safety;

        if let Some(kappa) = self.config.limiter {
            // Soederlind's smooth limiter keeps the step size differentiable in
            // the error, which avoids the on/off oscillation of hard clipping.
            raw = 1.0 + kappa * ((raw - 1.0) / kappa).atan();
        }

        let upper = if self.rejected_in_a_row > 0 {
            self.config.max_scale_after_reject
        } else {
            self.config.max_scale
        };
        let scale = raw.clamp(self.config.min_scale, upper);

        if accepted {
            self.errors[1] = self.errors[0];
            self.errors[0] = e0;
            self.ratios[1] = self.ratios[0];
            self.ratios[0] = if step_ratio.is_finite() && step_ratio > 0.0 {
                step_ratio
            } else {
                1.0
            };
            self.history = (self.history + 1).min(2);
            self.rejected_in_a_row = 0;
        } else {
            self.rejected_in_a_row += 1;
        }

        Decision { accepted, scale }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integral_controller_reproduces_the_classic_formula() {
        let mut c = Controller::new(ControllerConfig::preset(ControllerPreset::I));
        let d = c.control(0.5, 1.0, 4);
        assert!(d.accepted);
        let expected = 0.9 * 0.5f64.powf(-1.0 / 5.0);
        assert!((d.scale - expected).abs() < 1e-12);
    }

    #[test]
    fn rejection_does_not_grow_the_step() {
        let mut c = Controller::new(ControllerConfig::preset(ControllerPreset::Pi4020));
        let d = c.control(4.0, 1.0, 4);
        assert!(!d.accepted);
        assert!(d.scale < 1.0);
        assert_eq!(c.rejections(), 1);
    }

    #[test]
    fn filter_uses_history_only_once_it_exists() {
        let mut c = Controller::new(ControllerConfig::preset(ControllerPreset::PiHairer));
        let first = c.control(0.25, 1.0, 4);
        // The first accepted step falls back to integral control.
        assert!((first.scale - 0.9 * 0.25f64.powf(-0.2)).abs() < 1e-12);
        let second = c.control(0.25, first.scale, 4);
        let expected = 0.9 * 0.25f64.powf(-0.7 / 5.0) * 0.25f64.powf(0.4 / 5.0);
        assert!((second.scale - expected).abs() < 1e-12);
    }

    #[test]
    fn every_preset_resolves_by_name() {
        for preset in ControllerPreset::all() {
            assert_eq!(ControllerPreset::from_name(preset.name()), Some(*preset));
        }
    }

    #[test]
    fn limiter_bounds_the_growth() {
        let mut cfg = ControllerConfig::preset(ControllerPreset::I);
        cfg.limiter = Some(1.0);
        let mut c = Controller::new(cfg);
        let d = c.control(1e-8, 1.0, 4);
        assert!(d.scale < 3.0);
    }
}
