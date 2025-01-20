mod pc_saft_gly;
mod pc_saft_mix;
mod pc_saft_pure;
pub use pc_saft_gly::PcSaftGlyPure;
pub use pc_saft_mix::PcSaftMix;
pub use pc_saft_pure::PcSaftPure;
/// PcSaftError
use thiserror::Error;
#[derive(Debug, Error)]
pub enum PcSaftErr {
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
