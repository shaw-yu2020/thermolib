/// GiiTerm
pub struct GiiTerm {
    eta0: (f64, f64),
    eta1: (f64, f64),
    eta2: (f64, f64),
    eta3: (f64, f64),
    eta4: (f64, f64),
}
impl GiiTerm {
    pub fn t0d0(&mut self, eta: f64) -> f64 {
        self.eta0(eta)
    }
    pub fn t0d1(&mut self, eta: f64) -> f64 {
        eta * self.eta1(eta)
    }
    pub fn t0d2(&mut self, eta: f64) -> f64 {
        eta.powi(2) * self.eta2(eta)
    }
    pub fn t0d3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.eta3(eta)
    }
    pub fn t0d4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.eta4(eta)
    }
    pub fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 * self.eta1(eta)
    }
    pub fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 * (self.eta1(eta) + eta * self.eta2(eta))
    }
    pub fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 * eta * (2.0 * self.eta2(eta) + eta * self.eta3(eta))
    }
    pub fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 * eta.powi(2) * (3.0 * self.eta3(eta) + eta * self.eta4(eta))
    }
    pub fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        eta2 * self.eta1(eta) + eta1.powi(2) * self.eta2(eta)
    }
}
impl GiiTerm {
    pub fn lngii_t0d0(&mut self, eta: f64) -> f64 {
        self.eta0(eta).ln()
    }
    pub fn lngii_t0d1(&mut self, eta: f64) -> f64 {
        eta / self.eta0(eta) * self.eta1(eta)
    }
    pub fn lngii_t0d2(&mut self, eta: f64) -> f64 {
        (eta / self.eta0(eta)).powi(2) * (self.eta2(eta) * self.eta0(eta) - self.eta1(eta).powi(2))
    }
    pub fn lngii_t0d3(&mut self, eta: f64) -> f64 {
        (eta / self.eta0(eta)).powi(3)
            * (self.eta3(eta) * self.eta0(eta).powi(2)
                - 3.0 * self.eta2(eta) * self.eta1(eta) * self.eta0(eta)
                + 2.0 * self.eta1(eta).powi(3))
    }
    pub fn lngii_t0d4(&mut self, eta: f64) -> f64 {
        (eta / self.eta0(eta)).powi(4)
            * (self.eta4(eta) * self.eta0(eta).powi(3)
                - 4.0 * self.eta3(eta) * self.eta1(eta) * self.eta0(eta).powi(2)
                - 3.0 * self.eta2(eta).powi(2) * self.eta0(eta).powi(2)
                + 12.0 * self.eta2(eta) * self.eta1(eta).powi(2) * self.eta0(eta)
                - 6.0 * self.eta1(eta).powi(4))
    }
    pub fn lngii_t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 / self.eta0(eta) * self.eta1(eta)
    }
    pub fn lngii_t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 / self.eta0(eta).powi(2)
            * (eta * (self.eta2(eta) * self.eta0(eta) - self.eta1(eta).powi(2))
                + self.eta1(eta) * self.eta0(eta))
    }
    pub fn lngii_t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 * eta / self.eta0(eta).powi(3)
            * (eta
                * (self.eta3(eta) * self.eta0(eta).powi(2)
                    - 3.0 * self.eta2(eta) * self.eta1(eta) * self.eta0(eta)
                    + 2.0 * self.eta1(eta).powi(3))
                + 2.0
                    * (self.eta2(eta) * self.eta0(eta).powi(2)
                        - self.eta1(eta).powi(2) * self.eta0(eta)))
    }
    pub fn lngii_t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        eta1 * eta.powi(2) / self.eta0(eta).powi(4)
            * (eta
                * (self.eta4(eta) * self.eta0(eta).powi(3)
                    - 4.0 * self.eta3(eta) * self.eta1(eta) * self.eta0(eta).powi(2)
                    - 3.0 * self.eta2(eta).powi(2) * self.eta0(eta).powi(2)
                    + 12.0 * self.eta2(eta) * self.eta1(eta).powi(2) * self.eta0(eta)
                    - 6.0 * self.eta1(eta).powi(4))
                + 3.0
                    * (self.eta3(eta) * self.eta0(eta).powi(3)
                        - 3.0 * self.eta2(eta) * self.eta1(eta) * self.eta0(eta).powi(2)
                        + 2.0 * self.eta1(eta).powi(3) * self.eta0(eta)))
    }
    pub fn lngii_t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        self.eta0(eta).powi(2).recip()
            * (eta2 * self.eta1(eta) * self.eta0(eta)
                + eta1.powi(2) * (self.eta2(eta) * self.eta0(eta) - self.eta1(eta).powi(2)))
    }
}
impl GiiTerm {
    pub fn new() -> Self {
        Self {
            eta0: (0.0, 0.0),
            eta1: (0.0, 0.0),
            eta2: (0.0, 0.0),
            eta3: (0.0, 0.0),
            eta4: (0.0, 0.0),
        }
    }
    fn eta0(&mut self, eta: f64) -> f64 {
        if eta != self.eta0.0 {
            self.eta0 = (eta, (1.0 - 0.5 * eta) / (1.0 - eta).powi(3))
        }
        self.eta0.1
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        if eta != self.eta1.0 {
            self.eta1 = (eta, (2.5 - eta) / (1.0 - eta).powi(4))
        }
        self.eta1.1
    }
    fn eta2(&mut self, eta: f64) -> f64 {
        if eta != self.eta2.0 {
            self.eta2 = (eta, (9.0 - 3.0 * eta) / (1.0 - eta).powi(5))
        }
        self.eta2.1
    }
    fn eta3(&mut self, eta: f64) -> f64 {
        if eta != self.eta3.0 {
            self.eta3 = (eta, (42.0 - 12.0 * eta) / (1.0 - eta).powi(6))
        }
        self.eta3.1
    }
    fn eta4(&mut self, eta: f64) -> f64 {
        if eta != self.eta4.0 {
            self.eta4 = (eta, (240.0 - 60.0 * eta) / (1.0 - eta).powi(7))
        }
        self.eta4.1
    }
}
