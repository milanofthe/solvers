//! Additive Runge-Kutta methods: two tableaux on one set of abscissae.
//!
//! An IMEX method is a pair, and the pair is the method. Neither half means
//! anything alone: the explicit one integrates the part of the right hand side
//! that is cheap and non stiff, the implicit one the part that is stiff, and
//! the stages are shared, so the two are evaluated at the same points and the
//! coupling between them is what has to be right.
//!
//! Sharing the stages is not decoration, it is the condition that makes the
//! pair a method at all: the abscissae of the two tableaux have to agree, and a
//! file whose halves disagree on them is refused here rather than analysed.

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
        for i in 0..explicit.stages {
            let (left, right) = (explicit.c[i].value(), implicit.c[i].value());
            if (left - right).abs() > 1e-12 * left.abs().max(1.0) {
                return Err(TableauError(format!(
                    "the halves disagree on abscissa {}: {left} and {right}",
                    i + 1
                )));
            }
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
