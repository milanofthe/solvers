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
//! Reference: J. C. Butcher, "Numerical Methods for Ordinary Differential
//! Equations", 3rd ed., Wiley 2016, doi:10.1002/9781119121534

use crate::method::RkTableau;
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
pub fn elementary_weight(tableau: &RkTableau, weights: &[Coeff], tree: &Tree) -> Coeff {
    let psi = stage_weights(tableau, tree);
    let mut acc = Coeff::zero();
    for j in 0..tableau.stages {
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

/// Result of verifying a tableau against the order conditions.
#[derive(Clone, Debug, Serialize)]
pub struct OrderReport {
    /// Highest order whose conditions all hold.
    pub order: usize,
    /// Same for the embedded weights, if there are any.
    pub embedded_order: Option<usize>,
    /// Largest `q` with `sum_j a_ij c_j^(q-1) = c_i^q / q` for every stage.
    pub stage_order: usize,
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
    tableau: &RkTableau,
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
            let weight = elementary_weight(tableau, weights, &tree);
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
pub fn conditions_at(tableau: &RkTableau, weights: &[Coeff], order: usize) -> Vec<Condition> {
    trees_of_order(order)
        .into_iter()
        .map(|tree| {
            let weight = elementary_weight(tableau, weights, &tree);
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
pub fn error_constant(tableau: &RkTableau, weights: &[Coeff], order: usize) -> f64 {
    conditions_at(tableau, weights, order)
        .iter()
        .map(|condition| condition.residual * condition.residual)
        .sum::<f64>()
        .sqrt()
}

/// Verify a tableau. `max_order` bounds the search; ten is enough for every
/// published method and keeps the tree count manageable.
pub fn verify(tableau: &RkTableau, max_order: usize) -> OrderReport {
    let (order, exact, failing) = attained_order(tableau, &tableau.b, max_order);

    let embedded_order = tableau
        .b_embedded
        .as_ref()
        .map(|be| attained_order(tableau, be, max_order).0);

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
        consistent_abscissae: consistent,
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

    #[test]
    fn densities_of_the_first_trees_are_right() {
        assert_eq!(trees_of_order(1)[0].density().value(), 1.0);
        assert_eq!(trees_of_order(2)[0].density().value(), 2.0);
        let third: Vec<f64> = trees_of_order(3).iter().map(|t| t.density().value()).collect();
        // The two trees of order three have densities 3 and 6.
        assert!(third.contains(&3.0) && third.contains(&6.0));
    }
}
