//! A unified framework for ordinary differential equation solvers.
//!
//! The design principle is that a method is *data*. Butcher tableaux and
//! multistep coefficient patterns live in JSON files together with the
//! reference to the paper they come from; the code is generic over that data
//! and shared across method classes. What follows from it:
//!
//! * one Runge-Kutta stepper covering explicit, diagonally implicit and fully
//!   implicit tableaux, dispatching on the sparsity of `A`,
//! * one multistep engine that derives its variable step coefficients from the
//!   order conditions, so BDF, Adams-Bashforth, Adams-Moulton and friends are
//!   the same code with a different free coefficient pattern,
//! * one set of nonlinear solvers and one parametrized error controller shared
//!   by every implicit method,
//! * analyses (order conditions, stability regions, work-precision) that read
//!   the same data and therefore stay in sync with the methods.

pub mod analysis;
pub mod control;
pub mod linalg;
pub mod method;
pub mod nonlinear;
pub mod num;
pub mod ode;
pub mod problem;
pub mod problems;
pub mod simd;

#[cfg(feature = "embedded-methods")]
pub mod embedded {
    include!(concat!(env!("OUT_DIR"), "/embedded_methods.rs"));
}

pub use method::{Method, MethodClass, MethodLibrary};
pub use problem::{Problem, Stats};

/// The method library shipped with the crate.
#[cfg(feature = "embedded-methods")]
pub fn library() -> MethodLibrary {
    MethodLibrary::embedded().expect("the embedded method library must be valid")
}

/// A shared instance of the embedded library.
///
/// Multistep families that cannot start themselves name a start up method by
/// id, and this is where that lookup happens without rebuilding the library on
/// every solver construction.
#[cfg(feature = "embedded-methods")]
pub fn shared_library() -> &'static MethodLibrary {
    static LIBRARY: std::sync::OnceLock<MethodLibrary> = std::sync::OnceLock::new();
    LIBRARY.get_or_init(library)
}
