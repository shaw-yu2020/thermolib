/// GiiPure
pub struct GiiPure {
    neg_sum_xm1: f64,
    eta0: (f64, f64),
    eta1: (f64, f64),
    eta2: (f64, f64),
    eta3: (f64, f64),
    eta4: (f64, f64),
    lngii_eta1: (f64, f64),
    lngii_eta2: (f64, f64),
    lngii_eta3: (f64, f64),
    lngii_eta4: (f64, f64),
}
impl GiiPure {
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
impl GiiPure {
    pub fn lngii_mu_k<'a>(
        &mut self,
        eta: f64,
        rho_num: f64,
        m1_k: &'a [f64],
        eta_k: &'a [f64],
    ) -> impl Iterator<Item = f64> + use<'a> {
        let coef_m1 = -self.lngii_eta0(eta);
        let coef_eta = rho_num * self.neg_sum_xm1 * self.lngii_eta1(eta);
        m1_k.iter()
            .zip(eta_k)
            .map(move |(m1_k, eta_k)| coef_m1 * m1_k + coef_eta * eta_k)
    }
    pub fn lngii_t0d0(&mut self, eta: f64) -> f64 {
        self.neg_sum_xm1 * self.lngii_eta0(eta)
    }
    pub fn lngii_t0d1(&mut self, eta: f64) -> f64 {
        self.neg_sum_xm1 * eta * self.lngii_eta1(eta)
    }
    pub fn lngii_t0d2(&mut self, eta: f64) -> f64 {
        self.neg_sum_xm1 * eta.powi(2) * self.lngii_eta2(eta)
    }
    pub fn lngii_t0d3(&mut self, eta: f64) -> f64 {
        self.neg_sum_xm1 * eta.powi(3) * self.lngii_eta3(eta)
    }
    pub fn lngii_t0d4(&mut self, eta: f64) -> f64 {
        self.neg_sum_xm1 * eta.powi(4) * self.lngii_eta4(eta)
    }
    pub fn lngii_t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        self.neg_sum_xm1 * eta1 * self.lngii_eta1(eta)
    }
    pub fn lngii_t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        self.neg_sum_xm1 * eta1 * (self.lngii_eta1(eta) + eta * self.lngii_eta2(eta))
    }
    pub fn lngii_t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        self.neg_sum_xm1 * eta1 * eta * (2.0 * self.lngii_eta2(eta) + eta * self.lngii_eta3(eta))
    }
    pub fn lngii_t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        self.neg_sum_xm1
            * eta1
            * eta.powi(2)
            * (3.0 * self.lngii_eta3(eta) + eta * self.lngii_eta4(eta))
    }
    pub fn lngii_t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        self.neg_sum_xm1 * (eta2 * self.lngii_eta1(eta) + eta1.powi(2) * self.lngii_eta2(eta))
    }
}
impl GiiPure {
    pub fn new(sum_xm1: f64) -> Self {
        Self {
            neg_sum_xm1: -sum_xm1,
            eta0: (0.0, 0.0),
            eta1: (0.0, 0.0),
            eta2: (0.0, 0.0),
            eta3: (0.0, 0.0),
            eta4: (0.0, 0.0),
            lngii_eta1: (0.0, 0.0),
            lngii_eta2: (0.0, 0.0),
            lngii_eta3: (0.0, 0.0),
            lngii_eta4: (0.0, 0.0),
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
    #[inline]
    fn lngii_eta0(&mut self, eta: f64) -> f64 {
        self.eta0(eta).ln()
    }
    fn lngii_eta1(&mut self, eta: f64) -> f64 {
        if eta != self.lngii_eta1.0 {
            self.lngii_eta1 = (eta, (5.0 - 2.0 * eta) / ((eta - 2.0) * (eta - 1.0)))
        }
        self.lngii_eta1.1
    }
    fn lngii_eta2(&mut self, eta: f64) -> f64 {
        if eta != self.lngii_eta2.0 {
            self.lngii_eta2 = (
                eta,
                (2.0 * eta.powi(2) + 11.0 - 10.0 * eta) / ((eta - 2.0) * (eta - 1.0)).powi(2),
            )
        }
        self.lngii_eta2.1
    }
    fn lngii_eta3(&mut self, eta: f64) -> f64 {
        if eta != self.lngii_eta3.0 {
            self.lngii_eta3 = (
                eta,
                (30.0 * eta.powi(2) + 46.0 - (4.0 * eta.powi(3) + 66.0 * eta))
                    / ((eta - 2.0) * (eta - 1.0)).powi(3),
            )
        }
        self.lngii_eta3.1
    }
    fn lngii_eta4(&mut self, eta: f64) -> f64 {
        if eta != self.lngii_eta4.0 {
            self.lngii_eta4 = (
                eta,
                (12.0 * eta.powi(4) + 396.0 * eta.powi(2) + 282.0
                    - (120.0 * eta.powi(3) + 552.0 * eta))
                    / ((eta - 2.0) * (eta - 1.0)).powi(4),
            )
        }
        self.lngii_eta4.1
    }
}
