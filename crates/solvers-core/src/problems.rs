//! A catalogue of test problems.
//!
//! These are the problems the cost and convergence analyses run on. The set
//! follows the usual benchmark suite so that the numbers are comparable to the
//! published ones: non stiff orbit and population problems, the classical stiff
//! chemistry systems, and two problems with a known closed form so convergence
//! can be measured against an exact solution rather than a reference run.
//!
//! References
//! ----------
//! * E. Hairer, S. P. Noersett, G. Wanner, "Solving Ordinary Differential
//!   Equations I", 2nd ed., Springer 1993, doi:10.1007/978-3-540-78862-1
//! * E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II",
//!   2nd ed., Springer 1996, doi:10.1007/978-3-642-05221-7
//!
//! Where a problem splits the way an additive method wants, it says so. The
//! split is never invented: it is the one the term structure already has, the
//! stiff linear part on one side and what is left on the other, which is the
//! situation IMEX methods were designed for. A problem that states no split
//! gives everything to the implicit half, and an additive method run on it is
//! its implicit tableau.

use crate::linalg::Matrix;
use crate::problem::Problem;

/// A problem with the metadata a benchmark needs.
pub trait TestProblem: Problem {
    fn id(&self) -> &'static str;
    fn name(&self) -> &'static str;
    fn description(&self) -> &'static str;
    fn t_span(&self) -> (f64, f64);
    fn y0(&self) -> Vec<f64>;
    /// Rough stiffness classification, used to group results.
    fn is_stiff(&self) -> bool;
    /// Closed form solution, when one exists.
    fn exact(&self, _t: f64) -> Option<Vec<f64>> {
        None
    }
}

/// Every problem in the catalogue.
pub fn catalog() -> Vec<Box<dyn TestProblem>> {
    vec![
        Box::new(ExponentialDecay { lambda: 1.0 }),
        Box::new(NonlinearDecay::default()),
        Box::new(Oscillator { omega: 3.0 }),
        Box::new(Kaps { epsilon: 1e-3 }),
        Box::new(ProtheroRobinson { lambda: -1e4 }),
        Box::new(LotkaVolterra::default()),
        Box::new(Brusselator::default()),
        Box::new(VanDerPol { mu: 1.0 }),
        Box::new(VanDerPolStiff { mu: 1000.0 }),
        Box::new(Robertson),
        Box::new(Hires),
        Box::new(Orego),
        Box::new(Arenstorf),
    ]
}

pub fn get(id: &str) -> Option<Box<dyn TestProblem>> {
    catalog().into_iter().find(|p| p.id() == id)
}

// ---------------------------------------------------------------------------
// Problems with a closed form solution
// ---------------------------------------------------------------------------

/// `y' = -lambda * y`, the reference case for order verification.
pub struct ExponentialDecay {
    pub lambda: f64,
}

impl Problem for ExponentialDecay {
    fn dim(&self) -> usize {
        1
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -self.lambda * y[0];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, _y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = -self.lambda;
    }
}

impl TestProblem for ExponentialDecay {
    fn id(&self) -> &'static str {
        "exp_decay"
    }
    fn name(&self) -> &'static str {
        "Exponential decay"
    }
    fn description(&self) -> &'static str {
        "Scalar linear decay with a closed form solution, the cleanest way to measure the order of a method."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 5.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
    fn exact(&self, t: f64) -> Option<Vec<f64>> {
        Some(vec![(-self.lambda * t).exp()])
    }
}

/// `y1' = -y1 + y2^2`, `y2' = -y2`, solved by `(2e^-t - e^-2t, e^-t)`.
///
/// A linear problem only measures the order the stability function agrees with
/// the exponential to, which for several high order methods is higher than
/// their true order. The quadratic coupling here makes the elementary
/// differentials nonzero so the full set of tree conditions is exercised, while
/// both eigenvalues stay at minus one so even a method with a small stability
/// region has room to work.
pub struct NonlinearDecay {
    pub end: f64,
}

impl Default for NonlinearDecay {
    fn default() -> Self {
        NonlinearDecay { end: 5.0 }
    }
}

impl Problem for NonlinearDecay {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -y[0] + y[1] * y[1];
        dy[1] = -y[1];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = -1.0;
        j[(0, 1)] = 2.0 * y[1];
        j[(1, 0)] = 0.0;
        j[(1, 1)] = -1.0;
    }
    // Decay implicitly, coupling explicitly: the pattern an IMEX method exists
    // for, with the linear part on the side that can absorb it and the
    // nonlinearity on the side that is cheap.
    fn has_splitting(&self) -> bool {
        true
    }
    fn rhs_explicit(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1] * y[1];
        dy[1] = 0.0;
    }
    fn rhs_implicit(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -y[0];
        dy[1] = -y[1];
    }
}

impl TestProblem for NonlinearDecay {
    fn id(&self) -> &'static str {
        "nonlinear_decay"
    }
    fn name(&self) -> &'static str {
        "Nonlinear decay"
    }
    fn description(&self) -> &'static str {
        "Quadratically coupled decay with a closed form solution. This is the one to measure a method's true nonlinear order on."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, self.end)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 1.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
    fn exact(&self, t: f64) -> Option<Vec<f64>> {
        Some(vec![2.0 * (-t).exp() - (-2.0 * t).exp(), (-t).exp()])
    }
}

/// Harmonic oscillator written as a first order system.
pub struct Oscillator {
    pub omega: f64,
}

impl Problem for Oscillator {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = -self.omega * self.omega * y[0];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, _y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = 0.0;
        j[(0, 1)] = 1.0;
        j[(1, 0)] = -self.omega * self.omega;
        j[(1, 1)] = 0.0;
    }
}

impl TestProblem for Oscillator {
    fn id(&self) -> &'static str {
        "oscillator"
    }
    fn name(&self) -> &'static str {
        "Harmonic oscillator"
    }
    fn description(&self) -> &'static str {
        "Undamped linear oscillator. Purely imaginary eigenvalues, so it probes the imaginary axis of the stability region."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 20.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 0.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
    fn exact(&self, t: f64) -> Option<Vec<f64>> {
        Some(vec![
            (self.omega * t).cos(),
            -self.omega * (self.omega * t).sin(),
        ])
    }
}

/// Kaps' problem. Stiff but with the closed form `y = (exp(-2t), exp(-t))`,
/// which makes it the standard probe for order reduction on stiff problems.
pub struct Kaps {
    pub epsilon: f64,
}

impl Problem for Kaps {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let inv = 1.0 / self.epsilon;
        dy[0] = -(inv + 2.0) * y[0] + inv * y[1] * y[1];
        dy[1] = y[0] - y[1] - y[1] * y[1];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        let inv = 1.0 / self.epsilon;
        j[(0, 0)] = -(inv + 2.0);
        j[(0, 1)] = 2.0 * inv * y[1];
        j[(1, 0)] = 1.0;
        j[(1, 1)] = -1.0 - 2.0 * y[1];
    }
    // Everything carrying `1/epsilon` is the fast half, and it is the whole of
    // the first equation's stiffness. What is left is order one in time.
    fn has_splitting(&self) -> bool {
        true
    }
    fn rhs_explicit(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -2.0 * y[0];
        dy[1] = y[0] - y[1] - y[1] * y[1];
    }
    fn rhs_implicit(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let inv = 1.0 / self.epsilon;
        dy[0] = inv * (y[1] * y[1] - y[0]);
        dy[1] = 0.0;
    }
}

impl TestProblem for Kaps {
    fn id(&self) -> &'static str {
        "kaps"
    }
    fn name(&self) -> &'static str {
        "Kaps problem"
    }
    fn description(&self) -> &'static str {
        "Singularly perturbed system whose exact solution does not depend on the stiffness parameter at all, which is what makes order reduction show up here."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 1.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 1.0]
    }
    fn is_stiff(&self) -> bool {
        true
    }
    fn exact(&self, t: f64) -> Option<Vec<f64>> {
        Some(vec![(-2.0 * t).exp(), (-t).exp()])
    }
}

/// Prothero-Robinson: `y' = lambda (y - g(t)) + g'(t)` with `g(t) = cos(t)`.
/// The exact solution is `g` itself for consistent initial data.
pub struct ProtheroRobinson {
    pub lambda: f64,
}

impl Problem for ProtheroRobinson {
    fn dim(&self) -> usize {
        1
    }
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = self.lambda * (y[0] - t.cos()) - t.sin();
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, _y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = self.lambda;
    }
    // The stiffness is entirely in the first term and the second is a smooth
    // forcing, so the split is the way the problem is written.
    fn has_splitting(&self) -> bool {
        true
    }
    fn rhs_explicit(&self, t: f64, _y: &[f64], dy: &mut [f64]) {
        dy[0] = -t.sin();
    }
    fn rhs_implicit(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = self.lambda * (y[0] - t.cos());
    }
}

impl TestProblem for ProtheroRobinson {
    fn id(&self) -> &'static str {
        "prothero_robinson"
    }
    fn name(&self) -> &'static str {
        "Prothero-Robinson"
    }
    fn description(&self) -> &'static str {
        "The canonical stiff scalar test. If a method still converges at its classical order here, it is stiffly accurate in practice."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 5.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0]
    }
    fn is_stiff(&self) -> bool {
        true
    }
    fn exact(&self, t: f64) -> Option<Vec<f64>> {
        Some(vec![t.cos()])
    }
}

// ---------------------------------------------------------------------------
// Non stiff problems
// ---------------------------------------------------------------------------

pub struct LotkaVolterra {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl Default for LotkaVolterra {
    fn default() -> Self {
        LotkaVolterra {
            a: 1.5,
            b: 1.0,
            c: 3.0,
            d: 1.0,
        }
    }
}

impl Problem for LotkaVolterra {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = self.a * y[0] - self.b * y[0] * y[1];
        dy[1] = -self.c * y[1] + self.d * y[0] * y[1];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = self.a - self.b * y[1];
        j[(0, 1)] = -self.b * y[0];
        j[(1, 0)] = self.d * y[1];
        j[(1, 1)] = -self.c + self.d * y[0];
    }
}

impl TestProblem for LotkaVolterra {
    fn id(&self) -> &'static str {
        "lotka_volterra"
    }
    fn name(&self) -> &'static str {
        "Lotka-Volterra"
    }
    fn description(&self) -> &'static str {
        "Predator prey system with a closed orbit. Smooth, non stiff, and unforgiving about phase error over many periods."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 15.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 1.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
}

pub struct Brusselator {
    pub a: f64,
    pub b: f64,
}

impl Default for Brusselator {
    fn default() -> Self {
        Brusselator { a: 1.0, b: 3.0 }
    }
}

impl Problem for Brusselator {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = self.a + y[0] * y[0] * y[1] - (self.b + 1.0) * y[0];
        dy[1] = self.b * y[0] - y[0] * y[0] * y[1];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = 2.0 * y[0] * y[1] - (self.b + 1.0);
        j[(0, 1)] = y[0] * y[0];
        j[(1, 0)] = self.b - 2.0 * y[0] * y[1];
        j[(1, 1)] = -y[0] * y[0];
    }
}

impl TestProblem for Brusselator {
    fn id(&self) -> &'static str {
        "brusselator"
    }
    fn name(&self) -> &'static str {
        "Brusselator"
    }
    fn description(&self) -> &'static str {
        "Autocatalytic reaction settling onto a limit cycle. Mildly stiff during the fast transients of each cycle."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 20.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.5, 3.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
}

pub struct VanDerPol {
    pub mu: f64,
}

impl Problem for VanDerPol {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = self.mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = 0.0;
        j[(0, 1)] = 1.0;
        j[(1, 0)] = -2.0 * self.mu * y[0] * y[1] - 1.0;
        j[(1, 1)] = self.mu * (1.0 - y[0] * y[0]);
    }
}

impl TestProblem for VanDerPol {
    fn id(&self) -> &'static str {
        "van_der_pol"
    }
    fn name(&self) -> &'static str {
        "Van der Pol (mu = 1)"
    }
    fn description(&self) -> &'static str {
        "Relaxation oscillator in its non stiff regime."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 20.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![2.0, 0.0]
    }
    fn is_stiff(&self) -> bool {
        false
    }
}

/// The stiff Van der Pol oscillator, the standard stress test for step size
/// control: the solution alternates between slow drift and near discontinuous
/// switches.
pub struct VanDerPolStiff {
    pub mu: f64,
}

impl Problem for VanDerPolStiff {
    // y'' = mu (1 - y^2) y' - y, the form the benchmark literature uses. The
    // relaxation period is about 1.61 mu, so the interval below covers roughly
    // two cycles.
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = self.mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = 0.0;
        j[(0, 1)] = 1.0;
        j[(1, 0)] = -2.0 * self.mu * y[0] * y[1] - 1.0;
        j[(1, 1)] = self.mu * (1.0 - y[0] * y[0]);
    }
}

impl TestProblem for VanDerPolStiff {
    fn id(&self) -> &'static str {
        "van_der_pol_stiff"
    }
    fn name(&self) -> &'static str {
        "Van der Pol (mu = 1000)"
    }
    fn description(&self) -> &'static str {
        "Stiff relaxation oscillator. Alternates between smooth drift and near discontinuous switching, which punishes a sluggish controller."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 3000.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![2.0, 0.0]
    }
    fn is_stiff(&self) -> bool {
        true
    }
}

pub struct Arenstorf;

impl Problem for Arenstorf {
    fn dim(&self) -> usize {
        4
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        const MU: f64 = 0.012277471;
        const NU: f64 = 1.0 - MU;
        let (y1, y2, y1p, y2p) = (y[0], y[1], y[2], y[3]);
        let d1 = ((y1 + MU) * (y1 + MU) + y2 * y2).powf(1.5);
        let d2 = ((y1 - NU) * (y1 - NU) + y2 * y2).powf(1.5);
        dy[0] = y1p;
        dy[1] = y2p;
        dy[2] = y1 + 2.0 * y2p - NU * (y1 + MU) / d1 - MU * (y1 - NU) / d2;
        dy[3] = y2 - 2.0 * y1p - NU * y2 / d1 - MU * y2 / d2;
    }
}

impl TestProblem for Arenstorf {
    fn id(&self) -> &'static str {
        "arenstorf"
    }
    fn name(&self) -> &'static str {
        "Arenstorf orbit"
    }
    fn description(&self) -> &'static str {
        "Restricted three body problem on a periodic orbit. The near passes at the two bodies force the step size to vary by orders of magnitude."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 17.065216560157962)
    }
    fn y0(&self) -> Vec<f64> {
        vec![0.994, 0.0, 0.0, -2.00158510637908252240537862224]
    }
    fn is_stiff(&self) -> bool {
        false
    }
}

// ---------------------------------------------------------------------------
// Stiff chemistry
// ---------------------------------------------------------------------------

/// Robertson's reaction kinetics. Three species over eleven decades in time.
pub struct Robertson;

impl Problem for Robertson {
    fn dim(&self) -> usize {
        3
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
        dy[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
        dy[2] = 3.0e7 * y[1] * y[1];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j[(0, 0)] = -0.04;
        j[(0, 1)] = 1.0e4 * y[2];
        j[(0, 2)] = 1.0e4 * y[1];
        j[(1, 0)] = 0.04;
        j[(1, 1)] = -1.0e4 * y[2] - 6.0e7 * y[1];
        j[(1, 2)] = -1.0e4 * y[1];
        j[(2, 0)] = 0.0;
        j[(2, 1)] = 6.0e7 * y[1];
        j[(2, 2)] = 0.0;
    }
}

impl TestProblem for Robertson {
    fn id(&self) -> &'static str {
        "robertson"
    }
    fn name(&self) -> &'static str {
        "Robertson"
    }
    fn description(&self) -> &'static str {
        "Reaction kinetics spanning eleven decades in time. The benchmark every stiff solver is measured on."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 1.0e5)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 0.0, 0.0]
    }
    fn is_stiff(&self) -> bool {
        true
    }
}

/// HIRES, a plant photosynthesis model with eight species.
pub struct Hires;

impl Problem for Hires {
    fn dim(&self) -> usize {
        8
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
        dy[1] = 1.71 * y[0] - 8.75 * y[1];
        dy[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
        dy[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
        dy[4] = -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6];
        dy[5] = -280.0 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
        dy[6] = 280.0 * y[5] * y[7] - 1.81 * y[6];
        dy[7] = -280.0 * y[5] * y[7] + 1.81 * y[6];
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        j.fill(0.0);
        j[(0, 0)] = -1.71;
        j[(0, 1)] = 0.43;
        j[(0, 2)] = 8.32;
        j[(1, 0)] = 1.71;
        j[(1, 1)] = -8.75;
        j[(2, 2)] = -10.03;
        j[(2, 3)] = 0.43;
        j[(2, 4)] = 0.035;
        j[(3, 1)] = 8.32;
        j[(3, 2)] = 1.71;
        j[(3, 3)] = -1.12;
        j[(4, 4)] = -1.745;
        j[(4, 5)] = 0.43;
        j[(4, 6)] = 0.43;
        j[(5, 3)] = 0.69;
        j[(5, 4)] = 1.71;
        j[(5, 5)] = -280.0 * y[7] - 0.43;
        j[(5, 6)] = 0.69;
        j[(5, 7)] = -280.0 * y[5];
        j[(6, 5)] = 280.0 * y[7];
        j[(6, 6)] = -1.81;
        j[(6, 7)] = 280.0 * y[5];
        j[(7, 5)] = -280.0 * y[7];
        j[(7, 6)] = 1.81;
        j[(7, 7)] = -280.0 * y[5];
    }
}

impl TestProblem for Hires {
    fn id(&self) -> &'static str {
        "hires"
    }
    fn name(&self) -> &'static str {
        "HIRES"
    }
    fn description(&self) -> &'static str {
        "Eight species model of photosynthesis in plants. Moderately stiff with a long quiescent tail."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 321.8122)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057]
    }
    fn is_stiff(&self) -> bool {
        true
    }
}

/// The Oregonator, a three variable model of the Belousov-Zhabotinsky reaction.
pub struct Orego;

impl Problem for Orego {
    fn dim(&self) -> usize {
        3
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        const S: f64 = 77.27;
        const Q: f64 = 8.375e-6;
        const W: f64 = 0.161;
        dy[0] = S * (y[1] + y[0] * (1.0 - Q * y[0] - y[1]));
        dy[1] = (y[2] - (1.0 + y[0]) * y[1]) / S;
        dy[2] = W * (y[0] - y[2]);
    }
    fn has_analytic_jacobian(&self) -> bool {
        true
    }
    fn jacobian(&self, _t: f64, y: &[f64], j: &mut Matrix<f64>) {
        const S: f64 = 77.27;
        const Q: f64 = 8.375e-6;
        const W: f64 = 0.161;
        j[(0, 0)] = S * (1.0 - 2.0 * Q * y[0] - y[1]);
        j[(0, 1)] = S * (1.0 - y[0]);
        j[(0, 2)] = 0.0;
        j[(1, 0)] = -y[1] / S;
        j[(1, 1)] = -(1.0 + y[0]) / S;
        j[(1, 2)] = 1.0 / S;
        j[(2, 0)] = W;
        j[(2, 1)] = 0.0;
        j[(2, 2)] = -W;
    }
}

impl TestProblem for Orego {
    fn id(&self) -> &'static str {
        "orego"
    }
    fn name(&self) -> &'static str {
        "Oregonator"
    }
    fn description(&self) -> &'static str {
        "Belousov-Zhabotinsky reaction. Periodic with sharp spikes several decades apart in magnitude."
    }
    fn t_span(&self) -> (f64, f64) {
        (0.0, 360.0)
    }
    fn y0(&self) -> Vec<f64> {
        vec![1.0, 2.0, 3.0]
    }
    fn is_stiff(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn central_difference_jacobian<P: Problem + ?Sized>(
        problem: &P,
        t: f64,
        y: &[f64],
        j: &mut Matrix<f64>,
    ) {
        let n = problem.dim();
        let mut minus = vec![0.0; n];
        let mut plus = vec![0.0; n];
        let mut probe = y.to_vec();
        for column in 0..n {
            let delta = 1e-6 * y[column].abs().max(1.0);
            probe[column] = y[column] - delta;
            problem.rhs(t, &probe, &mut minus);
            probe[column] = y[column] + delta;
            problem.rhs(t, &probe, &mut plus);
            probe[column] = y[column];
            for row in 0..n {
                j[(row, column)] = (plus[row] - minus[row]) / (2.0 * delta);
            }
        }
    }

    #[test]
    fn analytic_jacobians_match_finite_differences() {
        for problem in catalog() {
            if !problem.has_analytic_jacobian() {
                continue;
            }
            let n = problem.dim();
            let y = problem.y0();
            let t = 0.1;
            let mut analytic = Matrix::zeros(n, n);
            let mut numeric = Matrix::zeros(n, n);
            problem.jacobian(t, &y, &mut analytic);
            // Central differences, because a one sided difference on a problem
            // with a Jacobian as steep as Robertson's is dominated by its own
            // truncation error rather than by any mistake in the analytic form.
            central_difference_jacobian(&*problem, t, &y, &mut numeric);
            for i in 0..n {
                let mut column_scale = 1.0f64;
                for j in 0..n {
                    column_scale = column_scale.max(analytic[(i, j)].abs());
                }
                for j in 0..n {
                    let scale = column_scale;
                    assert!(
                        (analytic[(i, j)] - numeric[(i, j)]).abs() < 1e-4 * scale,
                        "{} jacobian mismatch at ({i},{j}): {} vs {}",
                        problem.id(),
                        analytic[(i, j)],
                        numeric[(i, j)]
                    );
                }
            }
        }
    }

    #[test]
    fn exact_solutions_satisfy_the_equation() {
        // A closed form solution must have the derivative the right hand side
        // claims, checked by a central difference in time.
        for problem in catalog() {
            let Some(_) = problem.exact(0.0) else { continue };
            let t = 0.3;
            let d = 1e-6;
            let (a, b) = (problem.exact(t - d).unwrap(), problem.exact(t + d).unwrap());
            let mut dy = vec![0.0; problem.dim()];
            problem.rhs(t, &problem.exact(t).unwrap(), &mut dy);
            for i in 0..problem.dim() {
                let numeric = (b[i] - a[i]) / (2.0 * d);
                assert!(
                    (numeric - dy[i]).abs() < 1e-4 * dy[i].abs().max(1.0),
                    "{} exact solution inconsistent in component {i}: {numeric} vs {}",
                    problem.id(),
                    dy[i]
                );
            }
        }
    }
}
