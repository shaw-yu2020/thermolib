use super::Alpha0Dtau;
use serde::{Deserialize, Serialize};
///
/// 理想气体亥姆霍兹方程
///
#[derive(Serialize, Deserialize, Debug)]
struct IdealPolynomialTerm {
    flag: i32,
    a: f64,
    t: f64,
}
impl IdealPolynomialTerm {
    #[allow(non_snake_case)]
    fn calc(&self, dtau: &Alpha0Dtau, tau: f64, Tr: f64) -> f64 {
        match self.flag {
            1 => {
                self.a / tau.powf(self.t)
                    * match dtau {
                        Alpha0Dtau::D0 => -1.0,
                        Alpha0Dtau::D1 => self.t,
                        Alpha0Dtau::D2 => -self.t * (self.t + 1.0),
                    }
            }
            2 => {
                self.a * Tr.powf(self.t) / tau.powf(self.t)
                    * match dtau {
                        Alpha0Dtau::D0 => -1.0 / self.t / (self.t + 1.0),
                        Alpha0Dtau::D1 => 1.0 / (self.t + 1.0),
                        Alpha0Dtau::D2 => -1.0,
                    }
            }
            _ => {
                println!("no flag={} in ideal_polynomial_term\n", self.flag);
                0.0
            }
        }
    }
}
#[derive(Serialize, Deserialize, Debug)]
struct IdealPlankEinsteinTerm {
    flag: i32,
    v: f64,
    b: f64,
}
impl IdealPlankEinsteinTerm {
    #[allow(non_snake_case)]
    fn calc(&self, dtau: &Alpha0Dtau, tau: f64, Tr: f64) -> f64 {
        let b = match self.flag {
            1 => self.b,
            2 => self.b / Tr,
            _ => {
                println!("no flag={} in ideal_plank_einstein_term\n", self.flag);
                self.b
            }
        };
        let exp_bitau = (-b * tau).exp();
        self.v
            * match dtau {
                Alpha0Dtau::D0 => (1.0 - exp_bitau).ln(),
                Alpha0Dtau::D1 => exp_bitau * tau * b / (1.0 - exp_bitau),
                Alpha0Dtau::D2 => -exp_bitau * (tau * b / (1.0 - exp_bitau)).powi(2),
            }
    }
}
#[derive(Serialize, Deserialize, Debug)]
pub struct IdealHelmholtzEquation {
    a_1: f64,
    a_tau: f64,
    a_lntau: f64,
    poly_terms: Vec<IdealPolynomialTerm>,
    pe_terms: Vec<IdealPlankEinsteinTerm>,
}
impl IdealHelmholtzEquation {
    #[allow(non_snake_case)]
    pub fn calc(&self, dtau: Alpha0Dtau, tau: f64, Tr: f64) -> f64 {
        let mut alpha0d = 0.0;
        match dtau {
            Alpha0Dtau::D0 => {
                alpha0d += self.a_1 + self.a_tau * tau + self.a_lntau * tau.ln();
            }
            Alpha0Dtau::D1 => {
                alpha0d += self.a_tau * tau + self.a_lntau;
            }
            Alpha0Dtau::D2 => {
                alpha0d += -self.a_lntau;
            }
        };
        for term in self.poly_terms.iter() {
            alpha0d += term.calc(&dtau, tau, Tr);
        }
        for term in self.pe_terms.iter() {
            alpha0d += term.calc(&dtau, tau, Tr);
        }
        alpha0d
    }
}
