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
// Universal Model Constants
const A00: f64 = 0.9105631445;
const A01: f64 = 0.6361281449;
const A02: f64 = 2.6861347891;
const A03: f64 = -26.547362491;
const A04: f64 = 97.759208784;
const A05: f64 = -159.59154087;
const A06: f64 = 91.297774084;
const A10: f64 = -0.3084016918;
const A11: f64 = 0.1860531159;
const A12: f64 = -2.5030047259;
const A13: f64 = 21.419793629;
const A14: f64 = -65.255885330;
const A15: f64 = 83.318680481;
const A16: f64 = -33.746922930;
const A20: f64 = -0.0906148351;
const A21: f64 = 0.4527842806;
const A22: f64 = 0.5962700728;
const A23: f64 = -1.7241829131;
const A24: f64 = -4.1302112531;
const A25: f64 = 13.776631870;
const A26: f64 = -8.6728470368;
const B00: f64 = 0.7240946941;
const B01: f64 = 2.2382791861;
const B02: f64 = -4.0025849485;
const B03: f64 = -21.003576815;
const B04: f64 = 26.855641363;
const B05: f64 = 206.55133841;
const B06: f64 = -355.60235612;
const B10: f64 = -0.5755498075;
const B11: f64 = 0.6995095521;
const B12: f64 = 3.8925673390;
const B13: f64 = -17.215471648;
const B14: f64 = 192.67226447;
const B15: f64 = -161.82646165;
const B16: f64 = -165.20769346;
const B20: f64 = 0.0976883116;
const B21: f64 = -0.2557574982;
const B22: f64 = -9.1558561530;
const B23: f64 = 20.642075974;
const B24: f64 = -38.804430052;
const B25: f64 = 93.626774077;
const B26: f64 = -29.666905585;
