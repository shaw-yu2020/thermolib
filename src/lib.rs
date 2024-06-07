//!
//! thermolib
//! =========
//!
//! An open-source library for
//! the calculation of fluid properties.
//!
//!
/// Python wrappers
mod python;
/// Vdw EOS
mod vdw;
pub use vdw::Vdw;
