use super::AlphaDtauDdeltaPlus;
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
    fn calc(&self, dd: AlphaDtauDdeltaPlus, tau: f64, delta: f64) -> f64 {
        let pt = self.n * delta.powf(self.d) * tau.powf(self.t);
        match dd {
            AlphaDtauDdeltaPlus::D00 => pt,
            AlphaDtauDdeltaPlus::D01 => self.d * pt,
            AlphaDtauDdeltaPlus::D02 => (self.d - 1.0) * self.d * pt,
            AlphaDtauDdeltaPlus::D10 => self.t * pt,
            AlphaDtauDdeltaPlus::D20 => (self.t - 1.0) * self.t * pt,
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
    fn calc(&self, dd: AlphaDtauDdeltaPlus, tau: f64, delta: f64) -> f64 {
        let et = self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp();
        match dd {
            AlphaDtauDdeltaPlus::D00 => et,
            AlphaDtauDdeltaPlus::D01 => (self.d - self.l * delta.powf(self.l)) * et,
            AlphaDtauDdeltaPlus::D02 => {
                ((self.d - 1.0) * self.d
                    - (2.0 * self.d + self.l - 1.0) * self.l * delta.powf(self.l)
                    + self.l.powi(2) * delta.powf(self.l).powi(2))
                    * et
            }
            AlphaDtauDdeltaPlus::D10 => self.t * et,
            AlphaDtauDdeltaPlus::D20 => (self.t - 1.0) * self.t * et,
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
    fn calc(&self, dd: AlphaDtauDdeltaPlus, tau: f64, delta: f64) -> f64 {
        let gt = self.n
            * delta.powf(self.d)
            * tau.powf(self.t)
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2));
        match dd {
            AlphaDtauDdeltaPlus::D00 => gt,
            AlphaDtauDdeltaPlus::D01 => {
                (self.d - 2.0 * self.eta * delta * (delta - self.epsilon)) * gt
            }
            AlphaDtauDdeltaPlus::D02 => {
                ((self.d - 1.0) * self.d
                    - 2.0 * self.eta * delta.powi(2)
                    - 4.0 * self.d * self.eta * delta * (delta - self.epsilon)
                    + 4.0 * self.eta.powi(2) * delta.powi(2) * (delta - self.epsilon).powi(2))
                    * gt
            }
            AlphaDtauDdeltaPlus::D10 => (self.t - 2.0 * self.beta * tau * (tau - self.gamma)) * gt,
            AlphaDtauDdeltaPlus::D20 => {
                ((self.t - 1.0) * self.t
                    - 2.0 * self.beta * tau.powi(2)
                    - 4.0 * self.t * self.beta * tau * (tau - self.gamma)
                    + 4.0 * self.beta.powi(2) * tau.powi(2) * (tau - self.gamma).powi(2))
                    * gt
            }
        }
    }
}
