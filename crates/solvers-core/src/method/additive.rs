//! Additive Runge-Kutta methods: two tableaux on one set of abscissae.
//!
//! An IMEX method is a pair, and the pair is the method. Neither half means
//! anything alone: the explicit one integrates the part of the right hand side
//! that is cheap and non stiff, the implicit one the part that is stiff, and
//! the stages are shared, so the two are evaluated at the same points and the
//! coupling between them is what has to be right.
//!
//! The stages are shared, the abscissae need not be. Each half evaluates its
//! own right hand side at its own `c_i`, which is the row sums of its own
//! matrix, and several published pairs use two different sets: the strong
//! stability preserving ones put the explicit half on `[0, 1]` while the
//! implicit half sits elsewhere. On an autonomous problem it makes no
//! difference at all, and the order conditions never mention `c`. Whether the
//! two agree is therefore reported and not required.

use super::rk::{RkTableau, RkTableauFile, TableauError};
use crate::num::Field;
use serde::{Deserialize, Serialize};

/// The two tableaux as they are written in a method file.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct AdditiveFile {
    /// The tableau applied to the non stiff half. Strictly lower triangular.
    pub explicit: RkTableauFile,
    /// The tableau applied to the stiff half. Diagonally implicit in every
    /// published pair, because a fully implicit half would defeat the purpose.
    pub implicit: RkTableauFile,
}

/// A validated additive pair.
#[derive(Clone, Debug)]
pub struct AdditiveTableau {
    pub explicit: RkTableau,
    pub implicit: RkTableau,
}

impl AdditiveTableau {
    pub fn from_file(file: &AdditiveFile) -> Result<AdditiveTableau, TableauError> {
        let explicit = RkTableau::from_file(&file.explicit)?;
        let implicit = RkTableau::from_file(&file.implicit)?;
        if explicit.stages != implicit.stages {
            return Err(TableauError(format!(
                "the halves have {} and {} stages, and an additive pair shares its stages",
                explicit.stages, implicit.stages
            )));
        }
        if !explicit.is_explicit() {
            return Err(TableauError(
                "the explicit half is not explicit".to_string(),
            ));
        }
        Ok(AdditiveTableau { explicit, implicit })
    }

    pub fn stages(&self) -> usize {
        self.explicit.stages
    }

    /// Whether the two halves evaluate at the same points in time.
    pub fn shares_abscissae(&self) -> bool {
        (0..self.stages()).all(|i| {
            let (left, right) = (self.explicit.c[i].value(), self.implicit.c[i].value());
            (left - right).abs() <= 1e-12 * left.abs().max(1.0)
        })
    }

    /// An additive pair estimates its error from both halves, so it is adaptive
    /// only when both carry an embedded weight vector.
    pub fn has_embedded(&self) -> bool {
        self.explicit.b_embedded.is_some() && self.implicit.b_embedded.is_some()
    }

    /// Stages that cost an implicit solve. The first is free wherever the
    /// implicit half starts explicitly, which every published pair does.
    pub fn implicit_stages(&self) -> usize {
        (0..self.stages())
            .filter(|&i| !self.implicit.a[(i, i)].is_zero())
            .count()
    }
}
