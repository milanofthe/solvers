//! Order conditions from rooted trees.
//!
//! The order of a Runge-Kutta method is not taken from the method file, it is
//! derived from the tableau. For every rooted tree `t` up to some order the
//! elementary weight must match the reciprocal of the density,
//!
//! ```text
//! Phi(t) = sum_j b_j * psi_j(t) = 1 / gamma(t)
//! psi_j(t) = prod_over_children_u ( sum_k a_jk * psi_k(u) ),  psi_j(leaf) = 1
//! gamma(t) = |t| * prod_over_children gamma(u)
//! ```
//!
//! The arithmetic runs in `Coeff`, so for a tableau published as exact
//! fractions the conditions are checked as identities and not against a
//! tolerance. Methods with irrational coefficients fall back to a numeric
//! comparison, and the report says which of the two happened.
//!
//! A Rosenbrock method is checked against the same trees, with a different
//! recursion. Writing every quantity in the step as a B-series and reading off
//! the coefficient of each tree turns
//!
//! ```text
//! K_i = h f(y + sum_j alpha_ij K_j) + h f'(y) sum_j gamma_ij K_j
//! ```
//!
//! into, for a tree `t` with subtrees `t_1..t_m`,
//!
//! ```text
//! k_i(t) = prod_r ( sum_j alpha_ij k_j(t_r) )
//!          + [m = 1] sum_j gamma_ij k_j(t_1)
//! ```
//!
//! The second term exists only for a root with a single child because `f'(y)`
//! is one derivative applied to one argument, so it can only produce an
//! elementary differential whose root has one branch. Setting `gamma` to zero
//! recovers the Runge-Kutta recursion, which is the check that the two agree.
//!
//! Reference: J. C. Butcher, "Numerical Methods for Ordinary Differential
//! Equations", 3rd ed., Wiley 2016, doi:10.1002/9781119121534

use crate::method::{RkTableau, RosenbrockTableau};
use crate::num::{Coeff, Field};
use serde::Serialize;

/// A rooted tree, stored as its (canonically ordered) list of subtrees.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Tree {
    pub children: Vec<Tree>,
}

impl Tree {
    pub fn leaf() -> Tree {
        Tree {
            children: Vec::new(),
        }
    }

    /// Number of nodes, which is the order of the condition it generates.
    pub fn order(&self) -> usize {
        1 + self.children.iter().map(|c| c.order()).sum::<usize>()
    }

    /// Butcher's density `gamma(t)`.
    pub fn density(&self) -> Coeff {
        let mut acc = Coeff::from_i64(self.order() as i64);
        for child in &self.children {
            acc = acc * child.density();
        }
        acc
    }

    /// Symmetry `sigma(t)`, needed for the elementary differential count.
    pub fn symmetry(&self) -> Coeff {
        let mut sorted = self.children.clone();
        sorted.sort();
        let mut acc = Coeff::one();
        let mut index = 0;
        while index < sorted.len() {
            let mut count = 1;
            while index + count < sorted.len() && sorted[index + count] == sorted[index] {
                count += 1;
            }
            let mut factorial = Coeff::one();
            for f in 2..=count {
                factorial = factorial * Coeff::from_i64(f as i64);
            }
            acc = acc * factorial * sorted[index].symmetry().powi(count as i32);
            index += count;
        }
        acc
    }

    /// Bracket notation, e.g. `[[t],t]`.
    pub fn to_string_compact(&self) -> String {
        if self.children.is_empty() {
            return "t".to_string();
        }
        let inner: Vec<String> = self.children.iter().map(|c| c.to_string_compact()).collect();
        format!("[{}]", inner.join(","))
    }
}

/// All rooted trees with exactly `order` nodes, in a canonical order.
pub fn trees_of_order(order: usize) -> Vec<Tree> {
    if order == 0 {
        return Vec::new();
    }
    let mut by_order: Vec<Vec<Tree>> = vec![Vec::new(); order + 1];
    by_order[1] = vec![Tree::leaf()];
    for n in 2..=order {
        let mut pool: Vec<Tree> = Vec::new();
        for smaller in by_order.iter().take(n) {
            pool.extend(smaller.iter().cloned());
        }
        pool.sort_by_key(|t| (t.order(), t.to_string_compact()));
        let mut result = Vec::new();
        for children in multisets(&pool, 0, n - 1) {
            result.push(Tree { children });
        }
        by_order[n] = result;
    }
    by_order.remove(order)
}

/// All rooted trees with up to `order` nodes.
pub fn trees_up_to(order: usize) -> Vec<Tree> {
    (1..=order).flat_map(trees_of_order).collect()
}

/// Multisets of trees from `pool[start..]` whose orders sum to `remaining`.
/// Indices are chosen non decreasing, which enumerates every multiset once.
fn multisets(pool: &[Tree], start: usize, remaining: usize) -> Vec<Vec<Tree>> {
    if remaining == 0 {
        return vec![Vec::new()];
    }
    let mut out = Vec::new();
    for index in start..pool.len() {
        let order = pool[index].order();
        if order > remaining {
            continue;
        }
        for rest in multisets(pool, index, remaining - order) {
            let mut combination = Vec::with_capacity(rest.len() + 1);
            combination.push(pool[index].clone());
            combination.extend(rest);
            out.push(combination);
        }
    }
    out
}

/// How a method turns a tree into one weight per stage. The order conditions
/// are the same for every method that has such a rule; only the rule differs.
pub trait StageWeights {
    fn stages(&self) -> usize;
    fn weights(&self, tree: &Tree) -> Vec<Coeff>;
}

impl StageWeights for RkTableau {
    fn stages(&self) -> usize {
        self.stages
    }
    fn weights(&self, tree: &Tree) -> Vec<Coeff> {
        stage_weights(self, tree)
    }
}

impl StageWeights for RosenbrockTableau {
    fn stages(&self) -> usize {
        self.stages
    }
    fn weights(&self, tree: &Tree) -> Vec<Coeff> {
        rosenbrock_stage_weights(self, tree)
    }
}

/// Stage weights `k_j(t)` for a Rosenbrock method.
fn rosenbrock_stage_weights(tableau: &RosenbrockTableau, tree: &Tree) -> Vec<Coeff> {
    let s = tableau.stages;
    let children: Vec<Vec<Coeff>> = tree
        .children
        .iter()
        .map(|child| rosenbrock_stage_weights(tableau, child))
        .collect();

    let mut k = vec![Coeff::one(); s];
    for i in 0..s {
        for child in &children {
            let mut acc = Coeff::zero();
            for j in 0..s {
                if tableau.alpha[(i, j)].is_zero() {
                    continue;
                }
                acc = acc + tableau.alpha[(i, j)] * child[j];
            }
            k[i] = k[i] * acc;
        }
    }

    // The Jacobian term. It reaches back to the stage's own diagonal, which is
    // not circular: the value it needs there belongs to a strictly smaller tree
    // and has already been computed.
    if children.len() == 1 {
        for i in 0..s {
            let mut acc = Coeff::zero();
            for j in 0..s {
                if tableau.gamma[(i, j)].is_zero() {
                    continue;
                }
                acc = acc + tableau.gamma[(i, j)] * children[0][j];
            }
            k[i] = k[i] + acc;
        }
    }
    k
}

/// Stage weights `psi_j(t)` for one tree.
fn stage_weights(tableau: &RkTableau, tree: &Tree) -> Vec<Coeff> {
    let s = tableau.stages;
    let mut psi = vec![Coeff::one(); s];
    for child in &tree.children {
        let child_psi = stage_weights(tableau, child);
        for j in 0..s {
            let mut acc = Coeff::zero();
            for k in 0..s {
                if tableau.a[(j, k)].is_zero() {
                    continue;
                }
                acc = acc + tableau.a[(j, k)] * child_psi[k];
            }
            psi[j] = psi[j] * acc;
        }
    }
    psi
}

/// Elementary weight `Phi(t)` for a given set of weights.
pub fn elementary_weight(method: &dyn StageWeights, weights: &[Coeff], tree: &Tree) -> Coeff {
    let psi = method.weights(tree);
    let mut acc = Coeff::zero();
    for j in 0..method.stages() {
        acc = acc + weights[j] * psi[j];
    }
    acc
}

/// One order condition and how badly it is violated.
#[derive(Clone, Debug, Serialize)]
pub struct Condition {
    pub tree: String,
    pub order: usize,
    /// `Phi(t)` as computed from the tableau.
    pub weight: f64,
    /// `1 / gamma(t)`, the value it should have.
    pub target: f64,
    pub residual: f64,
    /// True when the comparison was made in exact arithmetic.
    pub exact: bool,
    pub satisfied: bool,
}

/// How far each of the three simplifying assumptions holds.
///
/// ```text
/// B(p): sum_i b_i c_i^(k-1) = 1/k                       k = 1..p
/// C(q): sum_j a_ij c_j^(k-1) = c_i^k / k       every i, k = 1..q
/// D(r): sum_i b_i c_i^(k-1) a_ij = b_j (1 - c_j^k)/k    every j, k = 1..r
/// ```
///
/// `B` says the weights integrate a polynomial exactly, `C` says each stage
/// does, and `D` is the same statement for the adjoint method. They are worth
/// having because of Butcher's theorem: a method satisfying all three has order
/// at least `min(p, q + r + 1, 2q + 2)`, with no reference to trees at all.
///
/// That is the only affordable route to the order of a high stage collocation
/// method. Gauss-Legendre on eight nodes has order sixteen, and the rooted trees
/// up to that order number in the hundreds of thousands; these are a few hundred
/// sums.
///
/// Reference: E. Hairer, S. P. Noersett, G. Wanner, "Solving Ordinary
/// Differential Equations I", 2nd ed., Springer 1993, Theorem II.7.4.
#[derive(Copy, Clone, Debug, Serialize)]
pub struct SimplifyingAssumptions {
    pub quadrature: usize,
    pub stage: usize,
    pub adjoint: usize,
}

impl SimplifyingAssumptions {
    /// The order Butcher's theorem certifies from these three alone.
    pub fn certified_order(&self) -> usize {
        self.quadrature
            .min(self.stage + self.adjoint + 1)
            .min(2 * self.stage + 2)
    }
}

/// Result of verifying a tableau against the order conditions.
#[derive(Clone, Debug, Serialize)]
pub struct OrderReport {
    /// Highest order established, by whichever route got furthest.
    pub order: usize,
    /// Same for the embedded weights, if there are any.
    pub embedded_order: Option<usize>,
    /// Largest `q` with `sum_j a_ij c_j^(q-1) = c_i^q / q` for every stage.
    pub stage_order: usize,
    /// How far the simplifying assumptions hold.
    pub assumptions: SimplifyingAssumptions,
    /// True when the order came from the rooted tree conditions rather than
    /// from Butcher's theorem, which is the case for every method whose order
    /// is low enough for the trees to be enumerated.
    pub from_trees: bool,
    /// Whether `c` equals the row sums of `A`.
    pub consistent_abscissae: bool,
    /// True when every condition could be checked in exact arithmetic.
    pub exact: bool,
    /// The conditions of the first order that fails, for diagnosis.
    pub failing: Vec<Condition>,
}

fn satisfied(residual: Coeff, magnitude: f64) -> bool {
    match residual {
        Coeff::Exact(r) => r.is_zero(),
        Coeff::Approx(v) => v.abs() <= 1e-12 * magnitude.max(1.0),
    }
}

/// Check the conditions for one set of weights and return the attained order.
fn attained_order(
    method: &dyn StageWeights,
    weights: &[Coeff],
    max_order: usize,
) -> (usize, bool, Vec<Condition>) {
    let mut order = 0;
    let mut exact = true;
    let mut failing = Vec::new();

    for n in 1..=max_order {
        let mut all_ok = true;
        let mut level = Vec::new();
        for tree in trees_of_order(n) {
            let weight = elementary_weight(method, weights, &tree);
            let target = Coeff::one() / tree.density();
            let residual = weight - target;
            let is_exact = weight.is_exact() && target.is_exact();
            exact &= is_exact;
            let ok = satisfied(residual, target.value().abs());
            all_ok &= ok;
            level.push(Condition {
                tree: tree.to_string_compact(),
                order: n,
                weight: weight.value(),
                target: target.value(),
                residual: residual.value(),
                exact: is_exact,
                satisfied: ok,
            });
        }
        if all_ok {
            order = n;
        } else {
            failing = level.into_iter().filter(|c| !c.satisfied).collect();
            break;
        }
    }
    (order, exact, failing)
}

/// Every condition at one order, satisfied or not.
///
/// At the order just above the one a method achieves, these residuals are the
/// leading error coefficients: which elementary differential the method is
/// least accurate on, and how much of the local error each one carries.
pub fn conditions_at(method: &dyn StageWeights, weights: &[Coeff], order: usize) -> Vec<Condition> {
    trees_of_order(order)
        .into_iter()
        .map(|tree| {
            let weight = elementary_weight(method, weights, &tree);
            let target = Coeff::one() / tree.density();
            let residual = weight - target;
            let is_exact = weight.is_exact() && target.is_exact();
            Condition {
                satisfied: satisfied(residual, target.value().abs()),
                tree: tree.to_string_compact(),
                order,
                weight: weight.value(),
                target: target.value(),
                residual: residual.value(),
                exact: is_exact,
            }
        })
        .collect()
}

/// Euclidean norm of the residuals at one order, the usual scalar measure of
/// how large a method's leading error term is.
pub fn error_constant(method: &dyn StageWeights, weights: &[Coeff], order: usize) -> f64 {
    conditions_at(method, weights, order)
        .iter()
        .map(|condition| condition.residual * condition.residual)
        .sum::<f64>()
        .sqrt()
}

/// How far `B`, `C` and `D` hold for one set of weights.
///
/// `C` does not involve the weights, so it is the same for a method and its
/// embedded partner; the other two are not.
pub fn simplifying_assumptions(
    tableau: &RkTableau,
    weights: &[Coeff],
    max_order: usize,
) -> SimplifyingAssumptions {
    let s = tableau.stages;
    let c = &tableau.c;

    let mut quadrature = 0;
    for k in 1..=max_order {
        let mut sum = Coeff::zero();
        for i in 0..s {
            sum = sum + weights[i] * c[i].powi(k as i32 - 1);
        }
        let target = Coeff::one() / Coeff::from_i64(k as i64);
        if !satisfied(sum - target, target.value().abs()) {
            break;
        }
        quadrature = k;
    }

    let mut stage = 0;
    for k in 1..=max_order {
        let ok = (0..s).all(|i| {
            let mut sum = Coeff::zero();
            for j in 0..s {
                sum = sum + tableau.a[(i, j)] * c[j].powi(k as i32 - 1);
            }
            let target = c[i].powi(k as i32) / Coeff::from_i64(k as i64);
            satisfied(sum - target, target.value().abs().max(1.0))
        });
        if !ok {
            break;
        }
        stage = k;
    }

    let mut adjoint = 0;
    for k in 1..=max_order {
        let ok = (0..s).all(|j| {
            let mut sum = Coeff::zero();
            for i in 0..s {
                sum = sum + weights[i] * c[i].powi(k as i32 - 1) * tableau.a[(i, j)];
            }
            let target = weights[j] * (Coeff::one() - c[j].powi(k as i32))
                / Coeff::from_i64(k as i64);
            satisfied(sum - target, target.value().abs().max(1.0))
        });
        if !ok {
            break;
        }
        adjoint = k;
    }

    SimplifyingAssumptions {
        quadrature,
        stage,
        adjoint,
    }
}

/// Verify a tableau.
///
/// `max_order` bounds the search through the rooted trees, which is what limits
/// it: the number of trees at order `n` grows faster than any method's stage
/// count. Where the trees run out, Butcher's theorem takes over. The report says
/// which of the two established the order, because they are different kinds of
/// evidence and a reader is entitled to know which one is on offer.
pub fn verify(tableau: &RkTableau, max_order: usize) -> OrderReport {
    let (tree_order, exact, failing) = attained_order(tableau, &tableau.b, max_order);

    // The assumptions are cheap, so they are checked well past where the trees
    // stop rather than only as far.
    let limit = 24;
    let assumptions = simplifying_assumptions(tableau, &tableau.b, limit);
    let certified = assumptions.certified_order();

    // A tree search that stopped short of its own limit found a condition that
    // genuinely fails, and that settles the order. One that ran to the limit
    // only means the search ended, and the assumptions may know more.
    let found_a_failure = tree_order < max_order;
    let order = if found_a_failure {
        tree_order
    } else {
        tree_order.max(certified)
    };

    let embedded_order = tableau.b_embedded.as_ref().map(|be| {
        let (embedded_tree, ..) = attained_order(tableau, be, max_order);
        if embedded_tree < max_order {
            embedded_tree
        } else {
            embedded_tree.max(simplifying_assumptions(tableau, be, limit).certified_order())
        }
    });

    // Abscissae consistency: c must be the row sums of A.
    let mut consistent = true;
    for i in 0..tableau.stages {
        let mut sum = Coeff::zero();
        for j in 0..tableau.stages {
            sum = sum + tableau.a[(i, j)];
        }
        if !satisfied(sum - tableau.c[i], 1.0) {
            consistent = false;
        }
    }

    // Stage order: the largest q for which every stage integrates polynomials
    // of degree q - 1 exactly. This is what governs order reduction on stiff
    // problems, so it is worth reporting next to the classical order.
    let mut stage_order = 0;
    for q in 1..=max_order {
        let mut ok = true;
        for i in 0..tableau.stages {
            let mut sum = Coeff::zero();
            for j in 0..tableau.stages {
                sum = sum + tableau.a[(i, j)] * tableau.c[j].powi(q as i32 - 1);
            }
            let target = tableau.c[i].powi(q as i32) / Coeff::from_i64(q as i64);
            if !satisfied(sum - target, target.value().abs()) {
                ok = false;
                break;
            }
        }
        if ok {
            stage_order = q;
        } else {
            break;
        }
    }

    OrderReport {
        order,
        embedded_order,
        stage_order,
        assumptions,
        from_trees: found_a_failure || certified <= tree_order,
        consistent_abscissae: consistent,
        exact,
        failing,
    }
}

/// Verify a Rosenbrock tableau.
///
/// Stage order is not reported: the notion belongs to a method whose stages are
/// approximations of the solution at their abscissae, and a Rosenbrock stage is
/// an increment, not a solution value. The abscissae are the row sums of
/// `alpha` by construction rather than by claim, so there is nothing to check
/// there either.
pub fn verify_rosenbrock(tableau: &RosenbrockTableau, max_order: usize) -> OrderReport {
    let (order, exact, failing) = attained_order(tableau, &tableau.b, max_order);
    let embedded_order = tableau
        .b_embedded
        .as_ref()
        .map(|be| attained_order(tableau, be, max_order).0);

    OrderReport {
        order,
        embedded_order,
        stage_order: 0,
        assumptions: SimplifyingAssumptions {
            quadrature: 0,
            stage: 0,
            adjoint: 0,
        },
        from_trees: true,
        consistent_abscissae: true,
        exact,
        failing,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rooted_tree_counts_match_the_known_sequence() {
        // 1, 1, 2, 4, 9, 20, 48, 115, 286 rooted trees, OEIS A000081.
        let expected = [1, 1, 2, 4, 9, 20, 48, 115];
        for (index, count) in expected.iter().enumerate() {
            assert_eq!(trees_of_order(index + 1).len(), *count, "order {}", index + 1);
        }
    }

    /// With no Jacobian coupling a Rosenbrock method is an explicit Runge-Kutta
    /// method, and the two recursions have to agree tree for tree. This is the
    /// check that the B-series derivation above did not go wrong.
    #[test]
    fn the_rosenbrock_recursion_reduces_to_the_runge_kutta_one() {
        use crate::method::coeff_serde::CoeffValue;
        use crate::method::{RkTableauFile, RosenbrockFile, RosenbrockTableau};

        let value = |n: i128, d: i128| CoeffValue(Coeff::rational(n, d));
        // Classical fourth order Runge-Kutta, written as a Rosenbrock method
        // with a Jacobian coupling that only carries the diagonal it needs to
        // be well posed at all.
        let file = RosenbrockFile {
            alpha: vec![
                vec![],
                vec![value(1, 2)],
                vec![value(0, 1), value(1, 2)],
                vec![value(0, 1), value(0, 1), value(1, 1)],
            ],
            gamma: vec![
                vec![value(1, 1)],
                vec![value(0, 1), value(1, 1)],
                vec![value(0, 1), value(0, 1), value(1, 1)],
                vec![value(0, 1), value(0, 1), value(0, 1), value(1, 1)],
            ],
            b: vec![value(1, 6), value(1, 3), value(1, 3), value(1, 6)],
            b_embedded: None,
        };
        let rosenbrock = RosenbrockTableau::from_file(&file).unwrap();

        // The same tableau as a Runge-Kutta method.
        let rk_file = RkTableauFile {
            a: vec![
                vec![],
                vec![value(1, 2)],
                vec![value(0, 1), value(1, 2)],
                vec![value(0, 1), value(0, 1), value(1, 1)],
            ],
            b: vec![value(1, 6), value(1, 3), value(1, 3), value(1, 6)],
            c: None,
            b_embedded: None,
            dense_output: None,
        };
        let rk = RkTableau::from_file(&rk_file).unwrap();

        for tree in trees_up_to(5) {
            let mut stripped = rosenbrock.clone();
            for i in 0..stripped.stages {
                for j in 0..stripped.stages {
                    stripped.gamma[(i, j)] = Coeff::zero();
                }
            }
            let a = rosenbrock_stage_weights(&stripped, &tree);
            let b = stage_weights(&rk, &tree);
            for (x, y) in a.iter().zip(&b) {
                assert_eq!(x.value(), y.value(), "tree {}", tree.to_string_compact());
            }
        }
    }

    /// One stage with `gamma = 1/2` is second order and no more, which is the
    /// smallest nontrivial statement the Rosenbrock conditions make.
    #[test]
    fn the_one_stage_rosenbrock_method_is_second_order_only_at_one_half() {
        use crate::method::coeff_serde::CoeffValue;
        use crate::method::{RosenbrockFile, RosenbrockTableau};

        for (numerator, denominator, expected) in [(1, 2, 2), (1, 1, 1)] {
            let file = RosenbrockFile {
                alpha: vec![vec![]],
                gamma: vec![vec![CoeffValue(Coeff::rational(numerator, denominator))]],
                b: vec![CoeffValue(Coeff::rational(1, 1))],
                b_embedded: None,
            };
            let tableau = RosenbrockTableau::from_file(&file).unwrap();
            assert_eq!(
                verify_rosenbrock(&tableau, 6).order,
                expected,
                "gamma = {numerator}/{denominator}"
            );
        }
    }

    #[test]
    fn densities_of_the_first_trees_are_right() {
        assert_eq!(trees_of_order(1)[0].density().value(), 1.0);
        assert_eq!(trees_of_order(2)[0].density().value(), 2.0);
        let third: Vec<f64> = trees_of_order(3).iter().map(|t| t.density().value()).collect();
        // The two trees of order three have densities 3 and 6.
        assert!(third.contains(&3.0) && third.contains(&6.0));
    }
}
