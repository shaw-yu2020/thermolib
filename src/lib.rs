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
/// Python wrappers
#[cfg(feature = "with_pyo3")]
mod python;
/// Vdw EOS
mod vdw;
pub use vdw::Vdw;
/// Rk EOS
mod rk;
pub use rk::Rk;
/// Srk EOS
mod srk;
pub use srk::Srk;
/// Pr EOS
mod pr;
pub use pr::Pr;
mod pr_mix;
pub use pr_mix::PrMix;
/// Helmholtz EOS
mod helmholtz;
pub use helmholtz::Helmholtz;
/// liquid metals
mod liquid_metal;
pub use liquid_metal::LiquidMetal;
/// PC-SAFT EOS
mod pc_saft;
pub use pc_saft::PcSaftGlyPure;
pub use pc_saft::PcSaftMix;
pub use pc_saft::PcSaftPure;
pub use pc_saft::PcSaftYglPure;
/// Ideal Gas
mod ideal_gas;
pub use ideal_gas::IdealGas;
/// Fundamental Constants
mod f64consts {
    pub use std::f64::consts::{FRAC_PI_2, FRAC_PI_6, PI, SQRT_2};
    pub const SQRT2ADD1: f64 = SQRT_2 + 1.0;
    pub const SQRT2SUB1: f64 = SQRT_2 - 1.0;
    pub const C: f64 = 299792458.0; // CODATA2018 (speed of light in vacuum) m s^-1
    pub const H: f64 = 6.62607015E-34; // CODATA2018 (Planck constant) J Hz^-1
    pub const K: f64 = 1.380649e-23; // CODATA2018 (Boltzmann constant) J K^-1
    pub const FRAC_H_K: f64 = H / K; // K / Hz
    pub const NA: f64 = 6.02214076E23; // CODATA2018 (Avogadro constant) mol^-1
    pub const FRAC_NA_1E30: f64 = NA / 1E30; // const FRAC_NA_1E30: f64 = 6.02214076E-7;
    pub const R: f64 = 8.314462618; // CODATA2018 (molar gas constant) J mol^-1 K^-1
}
