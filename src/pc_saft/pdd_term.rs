use crate::f64consts::K;
use std::f64::consts::PI;
pub struct PddTerm {
    temp: f64,
    rho_num: f64,
    epsilon: f64,
    a2_coef: f64,
    a3_coef: f64,
    j2a: J2aTerm,
    j2b: J2bTerm,
    j3c: J3cTerm,
}
impl PddTerm {
    pub fn new(m: f64, sigma3: f64, epsilon: f64, u: f64) -> Self {
        let u2 = u.powi(2) * 1E-19;
        Self {
            temp: 0.0,
            rho_num: 0.0,
            epsilon,
            a2_coef: -PI / sigma3 * (u2 / m / K).powi(2),
            a3_coef: -4.0 / 3.0 * PI.powi(2) / sigma3 * (u2 / m / K).powi(3),
            j2a: J2aTerm::new(m),
            j2b: J2bTerm::new(m),
            j3c: J3cTerm::new(m),
        }
    }
}
fn_polar!(PddTerm);
/// J2aTerm :: Universal Model Constants
const A00: f64 = 0.3043504;
const A01: f64 = -0.1358588;
const A02: f64 = 1.4493329;
const A03: f64 = 0.3556977;
const A04: f64 = -2.0653308;
const A10: f64 = 0.9534641;
const A11: f64 = -1.8396383;
const A12: f64 = 2.0131180;
const A13: f64 = -7.3724958;
const A14: f64 = 8.2374135;
const A20: f64 = -1.1610080;
const A21: f64 = 4.5258607;
const A22: f64 = 0.9751222;
const A23: f64 = -12.281038;
const A24: f64 = 5.9397575;
/// J2bTerm :: Universal Model Constants
const B00: f64 = 0.2187939;
const B01: f64 = -1.1896431;
const B02: f64 = 1.1626889;
const B10: f64 = -0.5873164;
const B11: f64 = 1.2489132;
const B12: f64 = -0.5085280;
const B20: f64 = 3.4869576;
const B21: f64 = -14.915974;
const B22: f64 = 15.372022;
/// J3cTerm :: Universal Model Constants
const C00: f64 = -0.0646774;
const C01: f64 = 0.1975882;
const C02: f64 = -0.8087562;
const C03: f64 = 0.6902849;
const C10: f64 = -0.9520876;
const C11: f64 = 2.9924258;
const C12: f64 = -2.3802636;
const C13: f64 = -0.2701261;
const C20: f64 = -0.6260979;
const C21: f64 = 1.2924686;
const C22: f64 = 1.6542783;
const C23: f64 = -3.4396744;
