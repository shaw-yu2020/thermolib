/// HsTerm
pub struct HsTerm {
    sum_xm: f64,
    eta0: (f64, f64),
    eta1: (f64, f64),
    eta2: (f64, f64),
    eta3: (f64, f64),
    eta4: (f64, f64),
}
impl HsTerm {
    pub fn t0d0(&mut self, eta: f64) -> f64 {
        self.sum_xm * self.eta0(eta)
    }
    pub fn t0d1(&mut self, eta: f64) -> f64 {
        self.sum_xm * eta * self.eta1(eta)
    }
    pub fn t0d2(&mut self, eta: f64) -> f64 {
        self.sum_xm * eta.powi(2) * self.eta2(eta)
    }
    pub fn t0d3(&mut self, eta: f64) -> f64 {
        self.sum_xm * eta.powi(3) * self.eta3(eta)
    }
    pub fn t0d4(&mut self, eta: f64) -> f64 {
        self.sum_xm * eta.powi(4) * self.eta4(eta)
    }
    pub fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        self.sum_xm * eta1 * self.eta1(eta)
    }
    pub fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        self.sum_xm * eta1 * (self.eta1(eta) + eta * self.eta2(eta))
    }
    pub fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        self.sum_xm * eta1 * eta * (2.0 * self.eta2(eta) + eta * self.eta3(eta))
    }
    pub fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        self.sum_xm * eta1 * eta.powi(2) * (3.0 * self.eta3(eta) + eta * self.eta4(eta))
    }
    pub fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        self.sum_xm * (eta2 * self.eta1(eta) + eta1.powi(2) * self.eta2(eta))
    }
}
impl HsTerm {
    pub fn new(sum_xm: f64) -> Self {
        Self {
            sum_xm,
            eta0: (0.0, 0.0),
            eta1: (0.0, 0.0),
            eta2: (0.0, 0.0),
            eta3: (0.0, 0.0),
            eta4: (0.0, 0.0),
        }
    }
    fn eta0(&mut self, eta: f64) -> f64 {
        if eta != self.eta0.0 {
            self.eta0 = (eta, eta * (4.0 - 3.0 * eta) / (1.0 - eta).powi(2))
        }
        self.eta0.1
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        if eta != self.eta1.0 {
            self.eta1 = (eta, (4.0 - 2.0 * eta) / (1.0 - eta).powi(3))
        }
        self.eta1.1
    }
    fn eta2(&mut self, eta: f64) -> f64 {
        if eta != self.eta2.0 {
            self.eta2 = (eta, (10.0 - 4.0 * eta) / (1.0 - eta).powi(4))
        }
        self.eta2.1
    }
    fn eta3(&mut self, eta: f64) -> f64 {
        if eta != self.eta3.0 {
            self.eta3 = (eta, (36.0 - 12.0 * eta) / (1.0 - eta).powi(5))
        }
        self.eta3.1
    }
    fn eta4(&mut self, eta: f64) -> f64 {
        if eta != self.eta4.0 {
            self.eta4 = (eta, (168.0 - 48.0 * eta) / (1.0 - eta).powi(6))
        }
        self.eta4.1
    }
}
