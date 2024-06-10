use serde::Deserialize;
#[derive(Deserialize)]
pub struct GkhTerm {
    n: f64,
    d: f64,
    t: f64,
    eta: f64,
    epsilon: f64,
    beta: f64,
    gamma: f64,
    b: f64,
}
impl GkhTerm {
    pub fn tau0delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2)
                + 1.0 / (self.beta * (tau - self.gamma).powi(2) + self.b))
                .exp()
    }
    pub fn tau0delta1(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2)
                + 1.0 / (self.beta * (tau - self.gamma).powi(2) + self.b))
                .exp()
            * (self.d - 2.0 * self.eta * (delta - self.epsilon) * delta)
    }
    pub fn tau0delta2(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2)
                + 1.0 / (self.beta * (tau - self.gamma).powi(2) + self.b))
                .exp()
            * (self.d * (self.d - 1.0)
                - 4.0 * self.eta * (delta - self.epsilon) * self.d * delta
                - (2.0 * self.eta - 4.0 * self.eta.powi(2) * (delta - self.epsilon).powi(2))
                    * delta.powi(2))
    }
    pub fn tau1delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2)
                + 1.0 / (self.beta * (tau - self.gamma).powi(2) + self.b))
                .exp()
            * (self.t
                - 2.0 * self.beta * (tau - self.gamma) * tau
                    / (self.beta * (tau - self.gamma).powi(2) + self.b).powi(2))
    }
    pub fn tau1delta1(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2)
                + 1.0 / (self.beta * (tau - self.gamma).powi(2) + self.b))
                .exp()
            * (self.t
                - 2.0 * self.beta * (tau - self.gamma) * tau
                    / (self.beta * (tau - self.gamma).powi(2) + self.b).powi(2))
            * (self.d - 2.0 * self.eta * (delta - self.epsilon) * delta)
    }
    pub fn tau2delta0(&self, tau: f64, delta: f64) -> f64 {
        (self.n * delta.powf(self.d) * tau.powf(self.t))
            * (-self.eta * (delta - self.epsilon).powi(2)
                + 1.0 / (self.beta * (tau - self.gamma).powi(2) + self.b))
                .exp()
            * ((self.t
                - 2.0 * self.beta * (tau - self.gamma) * tau
                    / (self.beta * (tau - self.gamma).powi(2) + self.b).powi(2))
            .powi(2)
                - self.t
                - 2.0 * self.beta * tau.powi(2)
                    / (self.beta * (tau - self.gamma).powi(2) + self.b).powi(2)
                + 8.0 * self.beta.powi(2) * (tau - self.gamma).powi(2) * tau.powi(2)
                    / (self.beta * (tau - self.gamma).powi(2) + self.b).powi(3))
    }
}
