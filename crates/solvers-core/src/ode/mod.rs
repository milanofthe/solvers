//! Time integration.
//!
//! The driver is method agnostic: it owns the step size, the error controller
//! and the output, and asks a `Stepper` to attempt single steps. The steppers
//! are where the method classes differ, and there are only two of them because
//! the tableau and the multistep coefficient pattern already carry everything
//! that distinguishes the individual methods.

pub mod decouple;
pub mod lmm;
mod newton_matrix;
pub mod rk;
pub mod rosenbrock;

pub use decouple::{DecoupledLinear, StageDecoupling};
pub use lmm::LmmStepper;
pub use newton_matrix::NewtonMatrix;
pub use rk::RkStepper;
use rosenbrock::RosenbrockStepper;

use crate::control::{Controller, ControllerConfig};
use crate::method::{Method, MethodKind};
use crate::nonlinear::NonlinearConfig;
use crate::problem::{Problem, Stats};
use crate::simd;
use serde::{Deserialize, Serialize};

/// Result of a single attempted step.
#[derive(Copy, Clone, Debug)]
pub struct StepOutcome {
    /// False when the step could not be computed at all, for instance because
    /// the stage equations did not converge.
    pub ok: bool,
    /// Error estimate already normalized by the tolerances, so `<= 1` passes.
    pub error: f64,
}

impl StepOutcome {
    pub fn failed() -> StepOutcome {
        StepOutcome {
            ok: false,
            error: f64::INFINITY,
        }
    }
}

/// A single step of some integration method.
pub trait Stepper<P: Problem + ?Sized> {
    /// Order the error estimate converges with, used by the controller.
    fn control_order(&self) -> usize;

    /// Current accepted state.
    fn state(&self) -> &[f64];

    /// Candidate state produced by the last `attempt`.
    fn proposed(&self) -> &[f64];

    /// Prepare for integration from `(t, y)`.
    fn start(&mut self, problem: &P, stats: &mut Stats, t: f64, y: &[f64]);

    /// Try a step of size `h` from `t`.
    fn attempt(&mut self, problem: &P, stats: &mut Stats, t: f64, h: f64) -> StepOutcome;

    /// Accept the candidate.
    fn commit(&mut self, t: f64, h: f64);

    /// Discard the candidate.
    fn reject(&mut self, h: f64);

    /// Interpolate inside the last accepted step, `theta` in `[0, 1]`.
    fn interpolate(&self, theta: f64, out: &mut [f64]) -> bool;

    /// Upper bound on how fast the step size may grow. Multistep methods lose
    /// accuracy when the history spacing changes abruptly.
    fn max_growth(&self) -> f64 {
        f64::INFINITY
    }

    /// Chance for the stepper to shorten the step the driver proposed.
    ///
    /// Used by the multistep engine during start up, where the formula is still
    /// running at reduced order and needs a smaller step to keep its local
    /// error at the level the full order method will produce.
    fn step_limit(&self, h: f64) -> f64 {
        h
    }

    /// Whether this stepper produces a usable error estimate.
    fn is_adaptive(&self) -> bool;
}

/// Integration settings.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Options {
    pub rtol: f64,
    pub atol: f64,
    /// Initial step size. `None` selects one from the problem.
    pub h0: Option<f64>,
    pub h_min: f64,
    pub h_max: f64,
    pub max_steps: u64,
    /// Use the error estimate. Fixed step integration when false.
    pub adaptive: bool,
    pub controller: ControllerConfig,
    pub nonlinear: NonlinearConfig,
    /// Steps a Jacobian may be reused for before it is refreshed.
    pub max_jacobian_age: u32,
    /// Diagonalize the stage coupling of a fully implicit method instead of
    /// factoring the whole `s * n` system. Off only for testing the two paths
    /// against each other.
    pub decouple_stages: bool,
    /// Give a method with no embedded pair an error estimate by step doubling,
    /// so it can run adaptively at roughly three times the cost per step.
    pub richardson: bool,
    /// Times the solution is reported at. `None` reports every accepted step.
    pub t_eval: Option<Vec<f64>>,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            rtol: 1e-6,
            atol: 1e-9,
            h0: None,
            h_min: 0.0,
            h_max: f64::INFINITY,
            max_steps: 1_000_000,
            adaptive: true,
            controller: ControllerConfig::default(),
            nonlinear: NonlinearConfig::default(),
            max_jacobian_age: 5,
            decouple_stages: true,
            richardson: true,
            t_eval: None,
        }
    }
}

impl Options {
    pub fn with_tolerances(rtol: f64, atol: f64) -> Options {
        Options {
            rtol,
            atol,
            ..Default::default()
        }
    }

    pub fn fixed_step(h: f64) -> Options {
        Options {
            h0: Some(h),
            adaptive: false,
            ..Default::default()
        }
    }
}

/// How the integration ended.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    Success,
    MaxStepsExceeded,
    StepSizeUnderflow,
    StepFailed,
}

/// Result of an integration run.
#[derive(Clone, Debug)]
pub struct Solution {
    pub t: Vec<f64>,
    pub y: Vec<Vec<f64>>,
    pub stats: Stats,
    pub status: Status,
    /// Step sizes of the accepted steps, useful for diagnosing controllers.
    pub steps: Vec<f64>,
}

impl Solution {
    pub fn last(&self) -> Option<&Vec<f64>> {
        self.y.last()
    }

    pub fn succeeded(&self) -> bool {
        self.status == Status::Success
    }
}

/// Build the stepper a method needs.
pub fn stepper_for<P: Problem + ?Sized>(
    method: &Method,
    dim: usize,
    options: &Options,
) -> Box<dyn Stepper<P>> {
    match &method.kind {
        MethodKind::RungeKutta(tableau) => Box::new(RkStepper::new(
            tableau,
            method.declared_order.unwrap_or(1) as usize,
            method.declared_embedded_order.map(|v| v as usize),
            dim,
            options,
        )),
        MethodKind::LinearMultistep(family) => {
            Box::new(LmmStepper::with_startup(family, dim, options, startup_for(family, dim, options)))
        }
        MethodKind::Rosenbrock(tableau) => Box::new(RosenbrockStepper::new(
            tableau,
            method.declared_order.unwrap_or(1) as usize,
            method.declared_embedded_order.map(|v| v as usize),
            dim,
            options,
        )),
    }
}

/// Build the one step method a multistep family needs to fill its history.
///
/// Only families that cannot ramp their own order down to one need this; the
/// rest start themselves. The id comes from the method file, so the choice of
/// starter is data like everything else.
fn startup_for(family: &crate::method::LmmFamily, dim: usize, options: &Options) -> Option<RkStepper> {
    let id = family.startup.as_ref()?;
    #[cfg(feature = "embedded-methods")]
    {
        let method = crate::shared_library().get(id)?;
        let tableau = method.tableau()?;
        return Some(RkStepper::new(
            tableau,
            method.declared_order.unwrap_or(1) as usize,
            method.declared_embedded_order.map(|v| v as usize),
            dim,
            options,
        ));
    }
    #[cfg(not(feature = "embedded-methods"))]
    {
        let _ = (id, dim, options);
        None
    }
}

/// Hairer's starting step size heuristic.
///
/// Two trial evaluations give a local estimate of the second derivative, which
/// fixes a step that keeps the leading error term near the tolerance.
///
/// Reference: Hairer, Noersett, Wanner, "Solving ODEs I", 2nd ed., II.4.
pub fn initial_step<P: Problem + ?Sized>(
    problem: &P,
    stats: &mut Stats,
    t: f64,
    y: &[f64],
    order: usize,
    direction: f64,
    options: &Options,
) -> f64 {
    let n = y.len();
    let mut f0 = vec![0.0; n];
    stats.rhs_evals += 1;
    problem.rhs(t, y, &mut f0);

    let mut scale = vec![0.0; n];
    simd::error_scale(options.atol, options.rtol, y, y, &mut scale);

    let d0 = simd::weighted_rms(y, &scale);
    let d1 = simd::weighted_rms(&f0, &scale);
    let h0 = if d0 < 1e-5 || d1 < 1e-5 {
        1e-6
    } else {
        0.01 * d0 / d1
    };

    let mut y1 = vec![0.0; n];
    for i in 0..n {
        y1[i] = y[i] + direction * h0 * f0[i];
    }
    let mut f1 = vec![0.0; n];
    stats.rhs_evals += 1;
    problem.rhs(t + direction * h0, &y1, &mut f1);

    let mut diff = vec![0.0; n];
    for i in 0..n {
        diff[i] = f1[i] - f0[i];
    }
    let d2 = simd::weighted_rms(&diff, &scale) / h0;

    let h1 = if d1.max(d2) <= 1e-15 {
        (h0 * 1e-3).max(1e-6)
    } else {
        (0.01 / d1.max(d2)).powf(1.0 / (order.max(1) + 1) as f64)
    };

    (100.0 * h0).min(h1).clamp(1e-14, options.h_max)
}

/// Integrate `problem` over `t_span` with `method`.
pub fn integrate<P: Problem + ?Sized>(
    method: &Method,
    problem: &P,
    t_span: (f64, f64),
    y0: &[f64],
    options: &Options,
) -> Solution {
    let mut stepper = stepper_for::<P>(method, y0.len(), options);
    integrate_with(&mut *stepper, problem, t_span, y0, options)
}

/// Integrate with an already constructed stepper.
pub fn integrate_with<P: Problem + ?Sized>(
    stepper: &mut dyn Stepper<P>,
    problem: &P,
    t_span: (f64, f64),
    y0: &[f64],
    options: &Options,
) -> Solution {
    let (t0, t_end) = t_span;
    let direction = if t_end >= t0 { 1.0 } else { -1.0 };
    let n = y0.len();

    let mut stats = Stats::default();
    let mut controller = Controller::new(options.controller);
    stepper.start(problem, &mut stats, t0, y0);

    let adaptive = options.adaptive && stepper.is_adaptive();
    let mut h = match options.h0 {
        Some(h0) => h0.abs(),
        None => initial_step(
            problem,
            &mut stats,
            t0,
            y0,
            stepper.control_order(),
            direction,
            options,
        ),
    };
    h = h.min(options.h_max).min((t_end - t0).abs());
    if h <= 0.0 {
        h = 1e-6;
    }

    // Output handling. Either every accepted step, or the requested grid.
    let mut eval_points: Vec<f64> = options.t_eval.clone().unwrap_or_default();
    if !eval_points.is_empty() {
        eval_points.sort_by(|a, b| {
            if direction > 0.0 {
                a.total_cmp(b)
            } else {
                b.total_cmp(a)
            }
        });
    }
    let dense = !eval_points.is_empty();
    let mut next_eval = 0usize;

    let mut out_t: Vec<f64> = Vec::new();
    let mut out_y: Vec<Vec<f64>> = Vec::new();
    let mut steps: Vec<f64> = Vec::new();
    if !dense {
        out_t.push(t0);
        out_y.push(y0.to_vec());
    } else {
        while next_eval < eval_points.len()
            && (eval_points[next_eval] - t0) * direction <= 0.0
        {
            out_t.push(eval_points[next_eval]);
            out_y.push(y0.to_vec());
            next_eval += 1;
        }
    }

    let mut t = t0;
    let mut h_previous = h;
    let mut status = Status::Success;
    let mut consecutive_failures = 0u32;

    while (t_end - t) * direction > 0.0 {
        if stats.steps >= options.max_steps {
            status = Status::MaxStepsExceeded;
            break;
        }
        // `h` is the step the controller owns. What is actually taken can be
        // shorter, because the stepper may ask for a shorter start up step and
        // because the last step is trimmed to land on the end of the interval.
        // Neither of those is an accuracy signal, so they must not feed back
        // into `h`.
        let remaining = (t_end - t).abs();
        let h_try = stepper.step_limit(h).min(remaining);
        if h_try <= 0.0 {
            break;
        }

        stats.steps += 1;
        let signed = direction * h_try;
        let outcome = stepper.attempt(problem, &mut stats, t, signed);

        let decision = if adaptive {
            controller.control(outcome.error, h_try / h_previous, stepper.control_order())
        } else {
            crate::control::Decision {
                accepted: outcome.ok,
                scale: 1.0,
            }
        };

        if outcome.ok && decision.accepted {
            consecutive_failures = 0;
            let t_next = t + signed;

            if dense {
                while next_eval < eval_points.len()
                    && (eval_points[next_eval] - t_next) * direction <= 0.0
                {
                    let target = eval_points[next_eval];
                    let theta = if signed != 0.0 {
                        ((target - t) / signed).clamp(0.0, 1.0)
                    } else {
                        1.0
                    };
                    let mut buffer = vec![0.0; n];
                    if !stepper.interpolate(theta, &mut buffer) {
                        let a = stepper.state();
                        let b = stepper.proposed();
                        for i in 0..n {
                            buffer[i] = a[i] + theta * (b[i] - a[i]);
                        }
                    }
                    out_t.push(target);
                    out_y.push(buffer);
                    next_eval += 1;
                }
            }

            stepper.commit(t, signed);
            stats.accepted += 1;
            steps.push(signed);
            t = t_next;
            h_previous = h_try;

            if !dense {
                out_t.push(t);
                out_y.push(stepper.state().to_vec());
            }
        } else {
            stats.rejected += 1;
            stepper.reject(signed);
            consecutive_failures += 1;
            if consecutive_failures > 50 {
                status = Status::StepFailed;
                break;
            }
        }

        if outcome.ok {
            let growth = decision.scale.min(stepper.max_growth());
            h = (h * growth).min(options.h_max);
        } else {
            // A failed step says nothing about accuracy, so shrink decisively
            // rather than asking the controller.
            h *= 0.5;
        }
        if h < options.h_min.max(1e-14) {
            status = Status::StepSizeUnderflow;
            break;
        }
    }

    Solution {
        t: out_t,
        y: out_y,
        stats,
        status,
        steps,
    }
}
