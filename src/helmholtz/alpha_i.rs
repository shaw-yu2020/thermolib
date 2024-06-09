use serde::Deserialize;
#[derive(Deserialize)]
pub struct IdealHelmholtz {
    a_1: f64,
    a_tau: f64,
    a_lntau: f64,
    poly_terms: Vec<PolynomialTerm>,
    pe_terms: Vec<PlankEinsteinTerm>,
}
#[allow(non_snake_case)]
impl IdealHelmholtz {
    pub fn tau0(&self, tau: f64, Tc: f64) -> f64 {
        (self.a_1 + self.a_tau * tau + self.a_lntau * tau.ln())
            + self
                .poly_terms
                .iter()
                .map(|term| term.tau0(tau, Tc))
                .sum::<f64>()
            + self
                .pe_terms
                .iter()
                .map(|term| term.tau0(tau, Tc))
                .sum::<f64>()
    }
    pub fn tau1(&self, tau: f64, Tc: f64) -> f64 {
        (self.a_tau * tau + self.a_lntau)
            + self
                .poly_terms
                .iter()
                .map(|term| term.tau1(tau, Tc))
                .sum::<f64>()
            + self
                .pe_terms
                .iter()
                .map(|term| term.tau1(tau, Tc))
                .sum::<f64>()
    }
    pub fn tau2(&self, tau: f64, Tc: f64) -> f64 {
        -self.a_lntau
            + self
                .poly_terms
                .iter()
                .map(|term| term.tau2(tau, Tc))
                .sum::<f64>()
            + self
                .pe_terms
                .iter()
                .map(|term| term.tau2(tau, Tc))
                .sum::<f64>()
    }
}
#[derive(Deserialize)]
struct PolynomialTerm {
    is_alpha: bool,
    a: f64,
    t: f64,
}
#[allow(non_snake_case)]
impl PolynomialTerm {
    fn tau0(&self, tau: f64, Tc: f64) -> f64 {
        if self.is_alpha {
            -self.a / tau.powf(self.t)
        } else {
            -self.a * Tc.powf(self.t) / self.t / (self.t + 1.0) / tau.powf(self.t)
        }
    }
    fn tau1(&self, tau: f64, Tc: f64) -> f64 {
        if self.is_alpha {
            self.a * self.t / tau.powf(self.t)
        } else {
            self.a * Tc.powf(self.t) / (self.t + 1.0) / tau.powf(self.t)
        }
    }
    fn tau2(&self, tau: f64, Tc: f64) -> f64 {
        if self.is_alpha {
            -self.a * self.t * (self.t + 1.0) / tau.powf(self.t)
        } else {
            -self.a * Tc.powf(self.t) / tau.powf(self.t)
        }
    }
}
#[derive(Deserialize)]
struct PlankEinsteinTerm {
    is_alpha: bool,
    v: f64,
    b: f64,
}
#[allow(non_snake_case)]
impl PlankEinsteinTerm {
    fn tau0(&self, tau: f64, Tc: f64) -> f64 {
        let bitau = if self.is_alpha { self.b } else { self.b / Tc } * tau;
        let exp_bitau = (-bitau).exp();
        self.v * (1.0 - exp_bitau).ln()
    }
    fn tau1(&self, tau: f64, Tc: f64) -> f64 {
        let bitau = if self.is_alpha { self.b } else { self.b / Tc } * tau;
        let exp_bitau = (-bitau).exp();
        self.v * exp_bitau * bitau / (1.0 - exp_bitau)
    }
    fn tau2(&self, tau: f64, Tc: f64) -> f64 {
        let bitau = if self.is_alpha { self.b } else { self.b / Tc } * tau;
        let exp_bitau = (-bitau).exp();
        -self.v * exp_bitau * (bitau / (1.0 - exp_bitau)).powi(2)
    }
}
