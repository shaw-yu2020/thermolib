use crate::f64consts::{NA, R};
const FRAC_NA_1E30: f64 = NA / 1E30;
const FRAC_RE30_NA: f64 = R / FRAC_NA_1E30;
#[macro_use]
mod macros;
mod assoc_pure;
mod disp_term;
mod gii_term;
mod hs_term;
mod pc_saft_gly_pure;
mod pc_saft_pure;
mod polar_term;
mod s_pc_saft_mix2;
use assoc_pure::{AssocGlyPure, AssocPure};
use disp_term::DispTerm;
use gii_term::GiiPure;
use hs_term::HsPure;
use polar_term::PolarTerm;
// public
pub use pc_saft_gly_pure::PcSaftGlyPure;
pub use pc_saft_pure::PcSaftPure;
pub use s_pc_saft_mix2::SPcSaftMix2;
/// PcSaftError
use thiserror::Error;
#[derive(Debug, Error)]
enum PcSaftErr {
    #[error("c_flash diverge")]
    NotConvForC,
    #[error("t_flash diverge")]
    NotConvForT,
    #[error("tp_flash diverge")]
    NotConvForTP,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
    #[error("property only in double phase")]
    OnlyInDoublePhase,
    #[error("tpz_flash diverge")]
    NotConvForTPZ,
    #[error("tx_flash diverge")]
    NotConvForTX,
    #[error("ty_flash diverge")]
    NotConvForTY,
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pc_saft_pure() {
        // PcSaftPure::SO2
        let (m, sigma, epsilon) = (2.8611, 2.6826, 205.35);
        let mut fluid = PcSaftPure::new_fluid(m, sigma, epsilon);
        _fn_test!(fluid);
        // PcSaftPure::H2O
        let mut fluid = PcSaftPure::new_fluid(1.0656, 3.0007, 366.51);
        fluid.set_2B_assoc_term(0.034868, 2500.7); // kappa_AB epsilon_AB
        _fn_test!(fluid);
        // PcSaftPure::CO2
        let mut fluid = PcSaftPure::new_fluid(1.5131, 3.1869, 163.33);
        fluid.set_QQ_polar_term(4.4); // |Q|(DA)
        _fn_test!(fluid);
        // PcSaftPure::Acetone
        let mut fluid = PcSaftPure::new_fluid(2.7447, 3.2742, 232.99);
        fluid.set_DD_polar_term(2.88); // |u|(D)
        _fn_test!(fluid);
        // PcSaftGlyPure::CH3OH
        let mut fluid = PcSaftGlyPure::new_fluid(1.5255, 3.23, 188.9);
        fluid.set_2B_assoc_term(0.035176, 2899.5, 1.0, 1.0, 1.0);
        _fn_test!(fluid);
    }
}
