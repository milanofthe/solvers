//! Order conditions for additive Runge-Kutta methods.
//!
//! An IMEX method integrates one right hand side split in two,
//!
//! ```text
//! y' = f_E(t, y) + f_I(t, y)
//! ```
//!
//! with a pair of tableaux on the same abscissae, the explicit one on the part
//! that is cheap and the implicit one on the part that is stiff:
//!
//! ```text
//! Y_i     = y_n + h sum_j ( aE_ij f_E(Y_j) + aI_ij f_I(Y_j) )
//! y_{n+1} = y_n + h sum_i ( bE_i  f_E(Y_i) + bI_i  f_I(Y_i) )
//! ```
//!
//! Its order is not the order of either half. Each half has to reach the order
//! on its own, and on top of that every way the two can be interleaved has to
//! be right, which is what the extra conditions are about: they are indexed by
//! rooted trees whose nodes each carry one of the two colours, so a tree of
//! order `n` now stands for as many conditions as it has colourings.
//!
//! The weight of a two coloured tree is built the same way as for one colour,
//! reading the colour of each node to decide which matrix applies to it:
//!
//! ```text
//! psi_i(leaf of colour m)  = 1
//! psi_i(t)                 = prod over children c of ( A^(colour of c) psi(c) )_i
//! Phi(t)                   = sum_i b^(colour of root)_i psi_i(t)
//! ```
//!
//! and the condition is `Phi(t) = 1/gamma(t)`, with the same density as the one
//! coloured tree underneath, since the density counts nodes and not colours.
//!
//! With both halves equal this collapses to the classical conditions, which is
//! the check the tests make.

use super::order::{satisfied_residual, Condition, Tree};
use crate::method::AdditiveTableau;
use crate::num::{Coeff, Field};

/// A rooted tree whose nodes carry one of two colours.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ColouredTree {
    /// 0 is the explicit half, 1 the implicit one.
    pub colour: u8,
    pub children: Vec<ColouredTree>,
}

impl ColouredTree {
    pub fn order(&self) -> usize {
        1 + self.children.iter().map(|c| c.order()).sum::<usize>()
    }

    /// The uncoloured tree underneath, which is what carries the density.
    pub fn shape(&self) -> Tree {
        Tree {
            children: self.children.iter().map(|c| c.shape()).collect(),
        }
    }

    /// A name that reads as the tree: `[..]` for the explicit colour and `{..}`
    /// for the implicit one.
    pub fn label(&self) -> String {
        let (open, close) = if self.colour == 0 { ('[', ']') } else { ('{', '}') };
        let inner: String = self.children.iter().map(|c| c.label()).collect();
        format!("{open}{inner}{close}")
    }

}

/// How deep the coupling conditions are searched.
///
/// The two coloured trees grow about four and a half times per order where the
/// uncoloured ones grow two and a half, so a search that is cheap for a single
/// tableau is not cheap here. Eight is well past the published additive
/// methods, which stop at five.
pub const ADDITIVE_SEARCH_LIMIT: usize = 8;

/// Every two coloured tree of every order up to `order`, indexed by order.
///
/// Built the way the uncoloured ones are: a tree of order `n` is a coloured
/// root over a multiset of smaller trees whose orders sum to `n - 1`. The pool
/// is kept sorted so that choosing indices in non decreasing order enumerates
/// each multiset exactly once.
pub fn coloured_trees_by_order(order: usize) -> Vec<Vec<ColouredTree>> {
    let mut by_order: Vec<Vec<ColouredTree>> = vec![Vec::new()];
    for n in 1..=order {
        let level = level_of_order(&by_order, n);
        by_order.push(level);
    }
    by_order
}

/// The coloured trees with `n` nodes, from the levels below them.
fn level_of_order(by_order: &[Vec<ColouredTree>], n: usize) -> Vec<ColouredTree> {
    let mut pool: Vec<ColouredTree> = Vec::new();
    for smaller in by_order.iter().take(n) {
        pool.extend(smaller.iter().cloned());
    }
    // By order first, so the scan below can stop at the first tree too large to
    // fit. The derived ordering compares the colour first, which interleaves the
    // orders and would make that scan drop combinations.
    pool.sort_by(|a, b| a.order().cmp(&b.order()).then_with(|| a.cmp(b)));
    let mut level = Vec::new();
    for children in multisets(&pool, 0, n - 1) {
        for colour in 0..2u8 {
            level.push(ColouredTree {
                colour,
                children: children.clone(),
            });
        }
    }
    level.sort();
    level
}

fn multisets(pool: &[ColouredTree], start: usize, remaining: usize) -> Vec<Vec<ColouredTree>> {
    if remaining == 0 {
        return vec![Vec::new()];
    }
    let mut out = Vec::new();
    for index in start..pool.len() {
        let order = pool[index].order();
        // The pool is sorted by order first, so the first tree too large to fit
        // ends the scan.
        if order > remaining {
            break;
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

/// `psi(t)`, the stage weights of one coloured tree.
fn stage_weights(method: &AdditiveTableau, tree: &ColouredTree) -> Vec<Coeff> {
    let s = method.stages();
    let mut psi = vec![Coeff::one(); s];
    for child in &tree.children {
        let child_psi = stage_weights(method, child);
        let a = if child.colour == 0 {
            &method.explicit.a
        } else {
            &method.implicit.a
        };
        for j in 0..s {
            let mut acc = Coeff::zero();
            for k in 0..s {
                if a[(j, k)].is_zero() {
                    continue;
                }
                acc = acc + a[(j, k)] * child_psi[k];
            }
            psi[j] = psi[j] * acc;
        }
    }
    psi
}

/// `Phi(t)` and the size of the sum it came out of.
fn elementary_weight(method: &AdditiveTableau, tree: &ColouredTree) -> (Coeff, f64) {
    let psi = stage_weights(method, tree);
    let b = if tree.colour == 0 {
        &method.explicit.b
    } else {
        &method.implicit.b
    };
    let mut acc = Coeff::zero();
    let mut scale = 0.0;
    for j in 0..method.stages() {
        let term = b[j] * psi[j];
        acc = acc + term;
        scale += term.value().abs();
    }
    (acc, scale)
}

/// What an additive method attains, and where it first fails.
#[derive(Clone, Debug)]
pub struct AdditiveReport {
    pub order: usize,
    /// The order each half reaches on its own, which is an upper bound for the
    /// pair and usually larger than it.
    pub explicit_order: usize,
    pub implicit_order: usize,
    /// Whether every coefficient involved was exact.
    pub exact: bool,
    /// What the embedded pair attains, where there is one.
    pub embedded_order: Option<usize>,
    /// The conditions that fail at the first order that is not attained.
    pub failing: Vec<Condition>,
    /// Whether the two halves share their abscissae, which they have to.
    pub consistent_abscissae: bool,
}

/// How far the coupling conditions hold, ignoring what fails and why.
fn attained(method: &AdditiveTableau, max_order: usize) -> usize {
    let mut levels: Vec<Vec<ColouredTree>> = vec![Vec::new()];
    let mut order = 0;
    for n in 1..=max_order {
        levels.push(level_of_order(&levels, n));
        let mut all_ok = true;
        for tree in &levels[n] {
            let (weight, scale) = elementary_weight(method, tree);
            let target = Coeff::one() / tree.shape().density();
            if !satisfied_residual(weight - target, target.value().abs().max(scale)) {
                all_ok = false;
                break;
            }
        }
        if !all_ok {
            break;
        }
        order = n;
    }
    order
}

/// Attained order of an additive pair, from the coupling conditions.
pub fn verify(method: &AdditiveTableau, max_order: usize) -> AdditiveReport {
    let mut levels: Vec<Vec<ColouredTree>> = vec![Vec::new()];
    let mut order = 0;
    let mut exact = true;
    let mut failing = Vec::new();

    for n in 1..=max_order {
        levels.push(level_of_order(&levels, n));
        let mut all_ok = true;
        let mut level = Vec::new();
        for tree in &levels[n] {
            let (weight, scale) = elementary_weight(method, tree);
            let target = Coeff::one() / tree.shape().density();
            let residual = weight - target;
            let is_exact = weight.is_exact() && target.is_exact();
            exact &= is_exact;
            if satisfied_residual(residual, target.value().abs().max(scale)) {
                continue;
            }
            all_ok = false;
            level.push(Condition {
                tree: tree.label(),
                order: n,
                weight: weight.value(),
                target: target.value(),
                residual: residual.value(),
                exact: is_exact,
                satisfied: false,
            });
        }
        if all_ok {
            order = n;
        } else {
            failing = level;
            break;
        }
    }

    // The embedded pair is the same tableaux under the other weights, so it is
    // measured by the same conditions.
    let embedded_order = match (&method.explicit.b_embedded, &method.implicit.b_embedded) {
        (Some(be), Some(bi)) => {
            let mut pair = method.clone();
            pair.explicit.b = be.clone();
            pair.implicit.b = bi.clone();
            Some(attained(&pair, max_order))
        }
        _ => None,
    };

    let half = |single: &crate::method::RkTableau| super::order::verify(single, max_order).order;
    let mut consistent = method.explicit.stages == method.implicit.stages;
    if consistent {
        for i in 0..method.stages() {
            let d = (method.explicit.c[i].value() - method.implicit.c[i].value()).abs();
            if d > 1e-12 * method.explicit.c[i].value().abs().max(1.0) {
                consistent = false;
            }
        }
    }

    AdditiveReport {
        order,
        embedded_order,
        explicit_order: half(&method.explicit),
        implicit_order: half(&method.implicit),
        exact,
        failing,
        consistent_abscissae: consistent,
    }
}

/// Which of the failing conditions are couplings rather than a failure of one
/// half on its own. Reported separately because they say different things: a
/// half that fails is the wrong tableau, a coupling that fails is a pair that
/// was never meant to go together.
pub fn coupling_failures(failing: &[Condition]) -> usize {
    failing
        .iter()
        .filter(|condition| {
            let mixed = condition.tree.contains('[') && condition.tree.contains('{');
            mixed
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// With both halves the same tableau the coupling conditions collapse onto
    /// the classical ones, so the pair has to attain exactly what the single
    /// method does.
    #[test]
    fn a_pair_of_identical_halves_is_the_method_itself() {
        for id in ["rk4", "rkdp54", "esdirk43", "radau_iia_5", "heun3"] {
            let library = crate::MethodLibrary::embedded().unwrap();
            let method = library.get(id).unwrap();
            let tableau = method.tableau().unwrap().clone();
            let single = super::super::order::verify(&tableau, 8).order;
            let pair = AdditiveTableau {
                explicit: tableau.clone(),
                implicit: tableau,
            };
            let report = verify(&pair, 8);
            assert_eq!(report.order, single, "{id}");
            assert!(report.consistent_abscissae, "{id}");
        }
    }

    #[test]
    fn the_colouring_counts_are_the_published_ones() {
        // Counted by hand from the shapes: one node in two colours; the chain
        // of two in four; at order three the chain (eight) and the root with
        // two leaves (six); at order four sixteen, twelve, sixteen and eight.
        let levels = coloured_trees_by_order(4);
        let counts: Vec<usize> = (1..=4).map(|n| levels[n].len()).collect();
        assert_eq!(counts, vec![2, 4, 14, 52]);
    }
}
