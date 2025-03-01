//!
//! thermolib
//! =========
//!
//! An open-source library for
//! the calculation of fluid properties.
//!
//!
/// algorithms
mod algorithms;
/// PC-SAFT EOS
mod pc_saft;
pub use pc_saft::PcSaftGlyPure;
/// Python wrappers
#[cfg(feature = "with_pyo3")]
mod python;
/// Fundamental Constants
mod f64consts {
    pub const NA: f64 = 6.02214076E23; // CODATA2018 (Avogadro constant) mol^-1
    pub const R: f64 = 8.314462618; // CODATA2018 (molar gas constant) J mol^-1 K^-1
}
