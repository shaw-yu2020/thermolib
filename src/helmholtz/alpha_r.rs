use super::GkhTerm;
use serde::{Deserialize, Serialize};
#[derive(Debug, Deserialize, Serialize)]
pub struct ResidualHelmholtz {
    poly_terms: Vec<PolynomialTerm>,
    exp_terms: Vec<ExponentialTerm>,
    gauss_terms: Vec<GaussianTerm>,
    #[serde(default)]
    gkh_terms: Vec<GkhTerm>,
}
impl ResidualHelmholtz {
    pub fn tau0delta0(&self, tau: f64, delta: f64) -> f64 {
        self.poly_terms
            .iter()
            .map(|term| term.tau0delta0(tau, delta))
            .sum::<f64>()
            + self
                .exp_terms
                .iter()
                .map(|term| term.tau0delta0(tau, delta))
                .sum::<f64>()
            + self
                .gauss_terms
                .iter()
                .map(|term| term.tau0delta0(tau, delta))
                .sum::<f64>()
            + self
                .gkh_terms
                .iter()
                .map(|term| term.tau0delta0(tau, delta))
                .sum::<f64>()
    }
    pub fn tau0delta1(&self, tau: f64, delta: f64) -> f64 {
        self.poly_terms
            .iter()
            .map(|term| term.tau0delta1(tau, delta))
            .sum::<f64>()
            + self
                .exp_terms
                .iter()
                .map(|term| term.tau0delta1(tau, delta))
                .sum::<f64>()
            + self
                .gauss_terms
                .iter()
                .map(|term| term.tau0delta1(tau, delta))
                .sum::<f64>()
            + self
                .gkh_terms
                .iter()
                .map(|term| term.tau0delta1(tau, delta))
                .sum::<f64>()
    }
    pub fn tau0delta2(&self, tau: f64, delta: f64) -> f64 {
        self.poly_terms
            .iter()
            .map(|term| term.tau0delta2(tau, delta))
            .sum::<f64>()
            + self
                .exp_terms
                .iter()
                .map(|term| term.tau0delta2(tau, delta))
                .sum::<f64>()
            + self
                .gauss_terms
                .iter()
                .map(|term| term.tau0delta2(tau, delta))
                .sum::<f64>()
            + self
                .gkh_terms
                .iter()
                .map(|term| term.tau0delta2(tau, delta))
                .sum::<f64>()
    }
    pub fn tau1delta0(&self, tau: f64, delta: f64) -> f64 {
        self.poly_terms
            .iter()
            .map(|term| term.tau1delta0(tau, delta))
            .sum::<f64>()
            + self
                .exp_terms
                .iter()
                .map(|term| term.tau1delta0(tau, delta))
                .sum::<f64>()
            + self
                .gauss_terms
                .iter()
                .map(|term| term.tau1delta0(tau, delta))
                .sum::<f64>()
            + self
                .gkh_terms
                .iter()
                .map(|term| term.tau1delta0(tau, delta))
                .sum::<f64>()
    }
    pub fn tau1delta1(&self, tau: f64, delta: f64) -> f64 {
        self.poly_terms
            .iter()
            .map(|term| term.tau1delta1(tau, delta))
            .sum::<f64>()
            + self
                .exp_terms
                .iter()
                .map(|term| term.tau1delta1(tau, delta))
                .sum::<f64>()
            + self
                .gauss_terms
                .iter()
                .map(|term| term.tau1delta1(tau, delta))
                .sum::<f64>()
            + self
                .gkh_terms
                .iter()
                .map(|term| term.tau1delta1(tau, delta))
                .sum::<f64>()
    }
    pub fn tau2delta0(&self, tau: f64, delta: f64) -> f64 {
        self.poly_terms
            .iter()
            .map(|term| term.tau2delta0(tau, delta))
            .sum::<f64>()
            + self
                .exp_terms
                .iter()
                .map(|term| term.tau2delta0(tau, delta))
                .sum::<f64>()
            + self
                .gauss_terms
                .iter()
                .map(|term| term.tau2delta0(tau, delta))
                .sum::<f64>()
            + self
                .gkh_terms
                .iter()
                .map(|term| term.tau2delta0(tau, delta))
                .sum::<f64>()
    }
}
#[derive(Debug, Deserialize, Serialize)]
struct PolynomialTerm {
    n: f64,
    d: f64,
    t: f64,
}
impl PolynomialTerm {
    fn tau0delta0(&self, tau: f64, delta: f64) -> f64 {
        self.n * delta.powf(self.d) * tau.powf(self.t)
    }
    fn tau0delta1(&self, tau: f64, delta: f64) -> f64 {
        self.n * delta.powf(self.d) * tau.powf(self.t) * self.d
    }
    fn tau0delta2(&self, tau: f64, delta: f64) -> f64 {
        self.n * delta.powf(self.d) * tau.powf(self.t) * self.d * (self.d - 1.0)
    }
    fn tau1delta0(&self, tau: f64, delta: f64) -> f64 {
        self.n * delta.powf(self.d) * tau.powf(self.t) * self.t
    }
    fn tau1delta1(&self, tau: f64, delta: f64) -> f64 {
        self.n * delta.powf(self.d) * tau.powf(self.t) * self.t * self.d
    }
    fn tau2delta0(&self, tau: f64, delta: f64) -> f64 {
        self.n * delta.powf(self.d) * tau.powf(self.t) * self.t * (self.t - 1.0)
    }
}
#[derive(Debug, Deserialize, Serialize)]
struct ExponentialTerm {
    n: f64,
    d: f64,
    t: f64,
    l: f64,
}
impl ExponentialTerm {
    fn tau0delta0(&self, tau: f64, delta: f64) -> f64 {
        self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp()
    }
    fn tau0delta1(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp())
            * (self.d - self.l * delta.powf(self.l))
    }
    fn tau0delta2(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp())
            * ((self.d - 1.0) * self.d
                - (2.0 * self.d + self.l - 1.0) * self.l * delta.powf(self.l)
                + (self.l * delta.powf(self.l)).powi(2))
    }
    fn tau1delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp()) * self.t
    }
    fn tau1delta1(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp())
            * self.t
            * (self.d - self.l * delta.powf(self.l))
    }
    fn tau2delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t) * (-delta.powf(self.l)).exp())
            * self.t
            * (self.t - 1.0)
    }
}
#[derive(Debug, Deserialize, Serialize)]
struct GaussianTerm {
    n: f64,
    d: f64,
    t: f64,
    eta: f64,
    epsilon: f64,
    beta: f64,
    gamma: f64,
}
impl GaussianTerm {
    fn tau0delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2))
                .exp()
    }
    fn tau0delta1(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2))
                .exp()
            * (self.d - 2.0 * self.eta * delta * (delta - self.epsilon))
    }
    fn tau0delta2(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2))
                .exp()
            * ((self.d - 1.0) * self.d
                - 2.0 * self.eta * delta.powi(2)
                - 4.0 * self.d * self.eta * delta * (delta - self.epsilon)
                + 4.0 * (self.eta * delta * (delta - self.epsilon)).powi(2))
    }
    fn tau1delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2))
                .exp()
            * (self.t - 2.0 * self.beta * tau * (tau - self.gamma))
    }
    fn tau1delta1(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2))
                .exp()
            * (self.t - 2.0 * self.beta * tau * (tau - self.gamma))
            * (self.d - 2.0 * self.eta * delta * (delta - self.epsilon))
    }
    fn tau2delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2) - self.beta * (tau - self.gamma).powi(2))
                .exp()
            * ((self.t - 1.0) * self.t
                - 2.0 * self.beta * tau.powi(2)
                - 4.0 * self.t * self.beta * tau * (tau - self.gamma)
                + 4.0 * (self.beta * tau * (tau - self.gamma)).powi(2))
    }
}
