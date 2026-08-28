//! Measures what the stage decoupling actually buys, on a problem large enough
//! for the linear algebra to dominate.

use solvers_core::linalg::Matrix;
use solvers_core::method::MethodLibrary;
use solvers_core::ode::{self, Options};
use solvers_core::problem::Problem;
use std::time::Instant;

/// Semi discretized reaction diffusion on a line, a dense enough Jacobian to
/// make the factorization the dominant cost.
struct Diffusion {
    n: usize,
}

impl Problem for Diffusion {
    fn dim(&self) -> usize {
        self.n
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let n = self.n;
        let k = (n * n) as f64;
        for i in 0..n {
            let left = if i == 0 { 0.0 } else { y[i - 1] };
            let right = if i + 1 == n { 0.0 } else { y[i + 1] };
            dy[i] = k * (left - 2.0 * y[i] + right) + y[i] * (1.0 - y[i]);
        }
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        let n = self.n;
        let k = (n * n) as f64;
        j.fill(0.0);
        for i in 0..n {
            j[(i, i)] = -2.0 * k + 1.0 - 2.0 * y[i];
            if i > 0 {
                j[(i, i - 1)] = k;
            }
            if i + 1 < n {
                j[(i, i + 1)] = k;
            }
        }
    }
}

#[test]
fn stage_decoupling_reduces_the_cost() {
    let library = MethodLibrary::embedded().unwrap();
    let problem = Diffusion { n: 80 };
    let y0: Vec<f64> = (0..problem.n)
        .map(|i| (std::f64::consts::PI * (i + 1) as f64 / (problem.n + 1) as f64).sin())
        .collect();

    for id in ["radau_iia_5", "gauss_legendre_6"] {
        let method = library.get(id).unwrap();
        let mut coupled = Options::with_tolerances(1e-6, 1e-8);
        coupled.decouple_stages = false;
        coupled.max_steps = 50_000;
        let mut decoupled = coupled.clone();
        decoupled.decouple_stages = true;

        let start = Instant::now();
        let a = ode::integrate(method, &problem, (0.0, 0.05), &y0, &coupled);
        let coupled_time = start.elapsed();

        let start = Instant::now();
        let b = ode::integrate(method, &problem, (0.0, 0.05), &y0, &decoupled);
        let decoupled_time = start.elapsed();

        assert!(a.succeeded() && b.succeeded(), "{id}: integration failed");
        let difference = a
            .last()
            .unwrap()
            .iter()
            .zip(b.last().unwrap())
            .fold(0.0f64, |acc, (x, y)| acc.max((x - y).abs()));

        println!(
            "{id:<18} n={} coupled {:>8.1?} ({} steps)  decoupled {:>8.1?} ({} steps)  speedup {:.2}x  difference {:.2e}",
            problem.n,
            coupled_time,
            a.stats.accepted,
            decoupled_time,
            b.stats.accepted,
            coupled_time.as_secs_f64() / decoupled_time.as_secs_f64(),
            difference
        );
    }
}
