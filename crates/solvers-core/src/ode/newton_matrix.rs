//! The iteration matrix shared by the implicit steppers.
//!
//! Stiff integrators spend most of their time in `I - h*gamma*J`. Keeping the
//! Jacobian and its factorization in one place lets a singly diagonal method
//! factor once per step instead of once per stage, and lets a method with an
//! unchanged step size skip the factorization entirely.

use crate::linalg::{Lu, Matrix};
use crate::problem::{Problem, Stats};

pub struct NewtonMatrix {
    /// Current approximation of df/dy.
    pub jacobian: Matrix<f64>,
    lu: Option<Lu<f64>>,
    /// Value of `h * gamma` the current factorization belongs to.
    factored_scale: f64,
    /// Number of steps the Jacobian has been reused for.
    age: u32,
    valid: bool,
}

impl NewtonMatrix {
    pub fn new(dim: usize) -> NewtonMatrix {
        NewtonMatrix {
            jacobian: Matrix::zeros(dim, dim),
            lu: None,
            factored_scale: f64::NAN,
            age: 0,
            valid: false,
        }
    }

    pub fn age(&self) -> u32 {
        self.age
    }

    /// Force a fresh Jacobian on the next request.
    pub fn invalidate(&mut self) {
        self.valid = false;
        self.lu = None;
    }

    /// Force a refactorization but keep the Jacobian.
    pub fn invalidate_factorization(&mut self) {
        self.lu = None;
        self.factored_scale = f64::NAN;
    }

    pub fn is_factored_for(&self, scale: f64) -> bool {
        self.lu.is_some() && self.factored_scale == scale
    }

    /// Make sure `I - scale * J` is factored, refreshing `J` when it is stale.
    pub fn prepare<P: Problem + ?Sized>(
        &mut self,
        problem: &P,
        stats: &mut Stats,
        t: f64,
        y: &[f64],
        scale: f64,
        max_age: u32,
    ) -> bool {
        let refresh = !self.valid || self.age >= max_age;
        if refresh {
            stats.jacobian_evals += 1;
            if !problem.has_analytic_jacobian() {
                stats.rhs_evals += problem.dim() as u64 + 1;
            }
            problem.jacobian(t, y, &mut self.jacobian);
            self.valid = true;
            self.age = 0;
            self.lu = None;
        }
        if self.is_factored_for(scale) {
            return true;
        }
        let n = self.jacobian.rows();
        let mut m = Matrix::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                m[(i, j)] = -scale * self.jacobian[(i, j)];
            }
            m[(i, i)] += 1.0;
        }
        let lu = Lu::factor(m);
        if lu.is_singular() {
            self.lu = None;
            return false;
        }
        stats.lu_decompositions += 1;
        self.factored_scale = scale;
        self.lu = Some(lu);
        true
    }

    /// Count one more step of reuse.
    pub fn advance_age(&mut self) {
        if self.valid {
            self.age += 1;
        }
    }

    pub fn solve(&self, rhs: &mut [f64]) -> bool {
        match &self.lu {
            Some(lu) => lu.solve_in_place(rhs),
            None => false,
        }
    }
}
