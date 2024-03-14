use super::Alpha0Dtau;
use super::AlphaDD;
use serde::{Deserialize, Serialize};
///
/// 理想气体亥姆霍兹方程
///
enum Dtau {
    D0, // 对无量纲温度的零阶导数
    D1, // 对无量纲温度的一阶导数
    D2, // 对无量纲温度的二阶导数
}
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
    fn calc(&self, dtau: Dtau, tau: f64, Tr: f64) -> f64 {
        let b = match self.flag {
            1 => self.b,
            2 => self.b / Tr,
            _ => {
                println!("no flag={} in ideal_plank_einstein_term\n", self.flag);
                self.b
            }
        };
        let exp_bitau = (-b * tau).exp();
        match dtau {
            Dtau::D0 => self.v * (1.0 - exp_bitau).ln(),
            Dtau::D1 => self.v * b * exp_bitau / (1.0 - exp_bitau),
            Dtau::D2 => -self.v * b.powi(2) * exp_bitau / (1.0 - exp_bitau).powi(2),
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
    pub fn calc(&self, dd: AlphaDD, tau: f64, delta: f64, Tr: f64) -> f64 {
        let mut alpha0dd = 0.0;
        match dd {
            AlphaDD::D00 => {
                alpha0dd += self.a_1 + self.a_tau * tau + self.a_lntau * tau.ln() + delta.ln();
                for term in self.poly_terms.iter() {
                    alpha0dd += term.calc(&Alpha0Dtau::D0, tau, Tr);
                }
                for term in self.pe_terms.iter() {
                    alpha0dd += term.calc(Dtau::D0, tau, Tr);
                }
            }
            AlphaDD::D01 => alpha0dd = 1.0 / delta,
            AlphaDD::D02 => alpha0dd = -1.0 / delta.powi(2),
            AlphaDD::D10 => {
                alpha0dd += self.a_tau + self.a_lntau / tau;
                for term in self.poly_terms.iter() {
                    alpha0dd += term.calc(&Alpha0Dtau::D1, tau, Tr);
                }
                for term in self.pe_terms.iter() {
                    alpha0dd += term.calc(Dtau::D1, tau, Tr);
                }
            }
            AlphaDD::D11 => alpha0dd = 0.0,
            AlphaDD::D20 => {
                alpha0dd += -self.a_lntau / tau.powi(2);
                for term in self.poly_terms.iter() {
                    alpha0dd += term.calc(&Alpha0Dtau::D2, tau, Tr);
                }
                for term in self.pe_terms.iter() {
                    alpha0dd += term.calc(Dtau::D2, tau, Tr);
                }
            }
        };
        alpha0dd
    }
}
