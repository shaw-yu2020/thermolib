use crate::f64consts::{NA, R};
const FRAC_NA_1E30: f64 = NA / 1E30;
const FRAC_RE30_NA: f64 = R / FRAC_NA_1E30;
#[macro_use]
mod macros;
mod disp_term;
mod gii_term;
mod hs_term;
mod pc_saft_gly;
use disp_term::DispTerm;
use gii_term::GiiPure;
use hs_term::HsPure;
pub use pc_saft_gly::PcSaftGlyPure;
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
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pc_saft_pure() {
        // PcSaftGlyPure::SO2
        let (m, sigma, epsilon) = (2.8611, 2.6826, 205.35);
        let mut fluid = PcSaftGlyPure::new_fluid(m, sigma, epsilon);
        _fn_test!(fluid);
        // PcSaftGlyPure::CH3OH
        let (m, sigma, epsilon) = (1.5255, 3.23, 188.9);
        let mut fluid = PcSaftGlyPure::new_fluid(m, sigma, epsilon);
        fluid.set_2B_assoc_type(0.035176, 2899.5);
        _fn_test!(fluid);
    }
}
