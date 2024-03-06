use super::AlphaDD;
use serde::{Deserialize, Serialize};
///
/// 剩余气体亥姆霍兹方程
///
#[derive(Serialize, Deserialize, Debug)]
struct ResidualPolynomialTerm {
    n: f64,
    d: f64,
    t: f64,
}
impl ResidualPolynomialTerm {
    fn calc(&self, dd: &AlphaDD, tau: f64, delta: f64) -> f64 {
        let pt = self.n * delta.powf(self.d) * tau.powf(self.t);
        match dd {
            AlphaDD::D00 => pt,
            AlphaDD::D01 => self.d * pt,
            AlphaDD::D02 => (self.d - 1.0) * self.d * pt,
            AlphaDD::D10 => self.t * pt,
            AlphaDD::D11 => self.d * self.t * pt,
            AlphaDD::D20 => (self.t - 1.0) * self.t * pt,
        }
    }
}
#[derive(Serialize, Deserialize, Debug)]
struct ResidualExponentialTerm {
    n: f64,
    d: f64,
    t: f64,
    l: f64,
}
impl ResidualExponentialTerm {
    fn calc(&self, dd: &AlphaDD, tau: f64, delta: f64) -> f64 {
        let et = self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp();
        match dd {
            AlphaDD::D00 => et,
            AlphaDD::D01 => (self.d - self.l * delta.powf(self.l)) * et,
            AlphaDD::D02 => {
                ((self.d - 1.0) * self.d
                    - (2.0 * self.d + self.l - 1.0) * self.l * delta.powf(self.l)
                    + self.l.powi(2) * delta.powf(self.l).powi(2))
                    * et
            }
            AlphaDD::D10 => self.t * et,
            AlphaDD::D11 => (self.d - self.l * delta.powf(self.l)) * self.t * et,
            AlphaDD::D20 => (self.t - 1.0) * self.t * et,
        }
    }
}
#[derive(Serialize, Deserialize, Debug)]
struct ResidualGaussianTerm {
    n: f64,
    d: f64,
    t: f64,
    eta: f64,
    epsilon: f64,
    beta: f64,
    gamma: f64,
}
impl ResidualGaussianTerm {
    fn calc(&self, dd: &AlphaDD, tau: f64, delta: f64) -> f64 {
        let gt = self.n
            * delta.powf(self.d)
            * tau.powf(self.t)
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2));
        match dd {
            AlphaDD::D00 => gt,
            AlphaDD::D01 => (self.d - 2.0 * self.eta * delta * (delta - self.epsilon)) * gt,
            AlphaDD::D02 => {
                ((self.d - 1.0) * self.d
                    - 2.0 * self.eta * delta.powi(2)
                    - 4.0 * self.d * self.eta * delta * (delta - self.epsilon)
                    + 4.0 * self.eta.powi(2) * delta.powi(2) * (delta - self.epsilon).powi(2))
                    * gt
            }
            AlphaDD::D10 => (self.t - 2.0 * self.beta * tau * (tau - self.gamma)) * gt,
            AlphaDD::D11 => {
                (self.d - 2.0 * self.eta * delta * (delta - self.epsilon))
                    * (self.t - 2.0 * self.beta * tau * (tau - self.gamma))
                    * gt
            }
            AlphaDD::D20 => {
                ((self.t - 1.0) * self.t
                    - 2.0 * self.beta * tau.powi(2)
                    - 4.0 * self.t * self.beta * tau * (tau - self.gamma)
                    + 4.0 * self.beta.powi(2) * tau.powi(2) * (tau - self.gamma).powi(2))
                    * gt
            }
        }
    }
}
#[derive(Serialize, Deserialize, Debug)]
pub struct ResidualHelmholtzEquation {
    poly_terms: Vec<ResidualPolynomialTerm>,
    exp_terms: Vec<ResidualExponentialTerm>,
    gauss_terms: Vec<ResidualGaussianTerm>,
}
impl ResidualHelmholtzEquation {
    pub fn calc(&self, dd: AlphaDD, tau: f64, delta: f64) -> f64 {
        let mut alphardd = 0.0;
        for term in self.poly_terms.iter() {
            alphardd += term.calc(&dd, tau, delta);
        }
        for term in self.exp_terms.iter() {
            alphardd += term.calc(&dd, tau, delta);
        }
        for term in self.gauss_terms.iter() {
            alphardd += term.calc(&dd, tau, delta);
        }
        alphardd
    }
}
