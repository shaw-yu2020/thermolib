use crate::f64consts::K;
use std::f64::consts::PI;
pub struct PqqTerm {
    temp: f64,
    rho_num: f64,
    epsilon: f64,
    a2_coef: f64,
    a3_coef: f64,
    j2a: J2aTerm,
    j2b: J2bTerm,
    j3c: J3cTerm,
}
impl PqqTerm {
    pub fn new(m: f64, sigma3: f64, epsilon: f64, q: f64) -> Self {
        let q2 = q.powi(2) * 1E-19;
        Self {
            temp: 0.0,
            rho_num: 0.0,
            epsilon,
            a2_coef: -(0.75 * q2 / m / K / sigma3).powi(2) * PI / sigma3.cbrt(),
            a3_coef: (0.75 * PI).powi(2) * (q2 / m / K / sigma3).powi(3),
            j2a: J2aTerm::new(m),
            j2b: J2bTerm::new(m),
            j3c: J3cTerm::new(m),
        }
    }
}
fn_polar!(PqqTerm);
/// J2aTerm :: Universal Model Constants
const A00: f64 = 1.2378308;
const A01: f64 = 2.4355031;
const A02: f64 = 1.6330905;
const A03: f64 = -1.6118152;
const A04: f64 = 6.9771185;
const A10: f64 = 1.2854109;
const A11: f64 = -11.465615;
const A12: f64 = 22.086893;
const A13: f64 = 7.4691383;
const A14: f64 = -17.197772;
const A20: f64 = 1.7942954;
const A21: f64 = 0.7695103;
const A22: f64 = 7.2647923;
const A23: f64 = 94.486699;
const A24: f64 = -77.148458;
/// J2bTerm :: Universal Model Constants
const B00: f64 = 0.4542718;
const B01: f64 = -4.5016264;
const B02: f64 = 3.5858868;
const B10: f64 = -0.8137340;
const B11: f64 = 10.064030;
const B12: f64 = -10.876631;
const B20: f64 = 6.8682675;
const B21: f64 = -5.1732238;
const B22: f64 = -17.240207;
/// J3cTerm :: Universal Model Constants
const C00: f64 = -0.5000437;
const C01: f64 = 6.5318692;
const C02: f64 = -16.014780;
const C03: f64 = 14.425970;
const C10: f64 = 2.0002094;
const C11: f64 = -6.7838658;
const C12: f64 = 20.383246;
const C13: f64 = -10.895984;
const C20: f64 = 3.1358271;
const C21: f64 = 7.2475888;
const C22: f64 = 3.0759478;
const C23: f64 = 0.0;
