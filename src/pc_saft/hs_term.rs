use std::f64::consts::FRAC_PI_6;
/// HsTerm
pub struct HsTerm {
    sum_xm: f64,
}
impl HsTerm {
    pub fn new(sum_xm: f64) -> Self {
        Self { sum_xm }
    }
    pub fn mu_k<'a>(
        &self,
        rho_num: f64,
        (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64),
        m_k: &'a [f64],
        zeta1_k: &'a [f64],
        zeta2_k: &'a [f64],
        zeta3_k: &'a [f64],
    ) -> impl Iterator<Item = f64> + use<'a> {
        let coef_m = -(1.0 - zeta3t0).ln();
        let coef_zeta1 = FRAC_PI_6.recip() * 3.0 * zeta2t0 / (1.0 - zeta3t0);
        let coef_zeta2 = FRAC_PI_6.recip()
            * 3.0
            * (zeta1t0 / (1.0 - zeta3t0)
                + zeta2t0.powi(2) / zeta3t0
                    * ((1.0 - zeta3t0).powi(2).recip() + (1.0 - zeta3t0).ln() / zeta3t0));
        let coef_zeta3 = FRAC_PI_6.recip()
            * (3.0 * zeta1t0 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                + 2.0 * zeta2t0.powi(3) / zeta3t0
                    * ((1.0 - zeta3t0).powi(3).recip() - (1.0 - zeta3t0).ln() / zeta3t0.powi(2))
                - zeta2t0.powi(3) * (2.0 - zeta3t0) / (zeta3t0 * (1.0 - zeta3t0)).powi(2))
            + rho_num * self.sum_xm / (1.0 - zeta3t0);
        m_k.iter().zip(zeta1_k).zip(zeta2_k).zip(zeta3_k).map(
            move |(((m_k, zeta1_k), zeta2_k), zeta3_k)| {
                coef_m * m_k + coef_zeta1 * zeta1_k + coef_zeta2 * zeta2_k + coef_zeta3 * zeta3_k
            },
        )
    }
    pub fn t0d0(&self, rho_num: f64, (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64)) -> f64 {
        (FRAC_PI_6 * rho_num).recip()
            * (3.0 * zeta1t0 * zeta2t0 / (1.0 - zeta3t0)
                + zeta2t0.powi(3) / zeta3t0
                    * ((1.0 - zeta3t0).powi(2).recip() + (1.0 - zeta3t0).ln() / zeta3t0))
            - self.sum_xm * (1.0 - zeta3t0).ln()
    }
    pub fn t0d1(&self, rho_num: f64, (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64)) -> f64 {
        (FRAC_PI_6 * rho_num).recip()
            * (3.0 * zeta1t0 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                + (zeta2t0 / (1.0 - zeta3t0)).powi(3) * (3.0 - zeta3t0))
            + self.sum_xm * zeta3t0 / (1.0 - zeta3t0)
    }
    pub fn t0d2(&self, rho_num: f64, (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64)) -> f64 {
        (FRAC_PI_6 * rho_num).recip()
            * (6.0 * zeta1t0 * zeta2t0 * zeta3t0 / (1.0 - zeta3t0).powi(3)
                + zeta2t0.powi(3) * (3.0 + 4.0 * zeta3t0 - zeta3t0.powi(2))
                    / (1.0 - zeta3t0).powi(4))
            + self.sum_xm * (zeta3t0 / (1.0 - zeta3t0)).powi(2)
    }
    pub fn t0d3(&self, rho_num: f64, (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64)) -> f64 {
        (FRAC_PI_6 * rho_num).recip()
            * (18.0 * zeta1t0 * zeta2t0 * zeta3t0.powi(2) / (1.0 - zeta3t0).powi(4)
                + 2.0 * zeta2t0.powi(3) * zeta3t0 * (8.0 + 5.0 * zeta3t0 - zeta3t0.powi(2))
                    / (1.0 - zeta3t0).powi(5))
            + self.sum_xm * 2.0 * (zeta3t0 / (1.0 - zeta3t0)).powi(3)
    }
    pub fn t1d0(
        &self,
        rho_num: f64,
        (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64),
        (zeta1t1, zeta2t1, zeta3t1): (f64, f64, f64),
    ) -> f64 {
        (FRAC_PI_6 * rho_num).recip()
            * (3.0
                * (zeta1t1 * zeta2t0 / (1.0 - zeta3t0)
                    + zeta2t1
                        * (zeta1t0 / (1.0 - zeta3t0)
                            + zeta2t0.powi(2) / zeta3t0
                                * ((1.0 - zeta3t0).powi(2).recip()
                                    + (1.0 - zeta3t0).ln() / zeta3t0)))
                + zeta3t1
                    * (3.0 * zeta1t0 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                        + 2.0 * zeta2t0.powi(3) / zeta3t0
                            * ((1.0 - zeta3t0).powi(3).recip()
                                - (1.0 - zeta3t0).ln() / zeta3t0.powi(2))
                        - zeta2t0.powi(3) / (zeta3t0 * (1.0 - zeta3t0)).powi(2) * (2.0 - zeta3t0)))
            + self.sum_xm * zeta3t1 / (1.0 - zeta3t0)
    }
    pub fn t1d1(
        &self,
        rho_num: f64,
        (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64),
        (zeta1t1, zeta2t1, zeta3t1): (f64, f64, f64),
    ) -> f64 {
        (FRAC_PI_6 * rho_num).recip()
            * (3.0
                * (zeta1t1 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                    + zeta2t1
                        * (zeta1t0 / (1.0 - zeta3t0).powi(2)
                            + zeta2t0.powi(2) * (3.0 - zeta3t0) / (1.0 - zeta3t0).powi(3)))
                + zeta3t1
                    * (6.0 * zeta1t0 * zeta2t0 / (1.0 - zeta3t0).powi(3)
                        + zeta2t0.powi(3) * (8.0 - 2.0 * zeta3t0) / (1.0 - zeta3t0).powi(4)))
            + self.sum_xm * zeta3t1 / (1.0 - zeta3t0).powi(2)
    }
    pub fn t2d0(
        &self,
        rho_num: f64,
        (zeta1t0, zeta2t0, zeta3t0): (f64, f64, f64),
        (zeta1t1, zeta2t1, zeta3t1): (f64, f64, f64),
        (zeta1t2, zeta2t2, zeta3t2): (f64, f64, f64),
    ) -> f64 {
        (FRAC_PI_6 * rho_num).recip()
            * (3.0
                * (zeta1t2 * zeta2t0 / (1.0 - zeta3t0)
                    + zeta2t2
                        * (zeta1t0 / (1.0 - zeta3t0)
                            + zeta2t0.powi(2) / zeta3t0
                                * ((1.0 - zeta3t0).powi(2).recip()
                                    + (1.0 - zeta3t0).ln() / zeta3t0)))
                + zeta3t2
                    * (3.0 * zeta1t0 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                        + 2.0 * zeta2t0.powi(3) / zeta3t0
                            * ((1.0 - zeta3t0).powi(3).recip()
                                - (1.0 - zeta3t0).ln() / zeta3t0.powi(2))
                        - zeta2t0.powi(3) / (zeta3t0 * (1.0 - zeta3t0)).powi(2) * (2.0 - zeta3t0))
                + 6.0
                    * (zeta1t1 * zeta2t1 / (1.0 - zeta3t0)
                        + zeta1t1 * zeta3t1 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                        + zeta2t1.powi(2) * zeta2t0 / zeta3t0
                            * ((1.0 - zeta3t0).powi(2).recip() + (1.0 - zeta3t0).ln() / zeta3t0)
                        + zeta2t1
                            * zeta3t1
                            * (zeta1t0 / (1.0 - zeta3t0).powi(2)
                                + 2.0 * zeta2t0.powi(2) / zeta3t0
                                    * ((1.0 - zeta3t0).powi(3).recip()
                                        - (1.0 - zeta3t0).ln() / zeta3t0.powi(2))
                                - (zeta2t0 / zeta3t0 / (1.0 - zeta3t0)).powi(2) * (2.0 - zeta3t0)))
                + zeta3t1.powi(2)
                    * (6.0 * zeta1t0 * zeta2t0 / (1.0 - zeta3t0).powi(3)
                        + 6.0 * zeta2t0.powi(3) / zeta3t0 / (1.0 - zeta3t0).powi(4)
                        + zeta2t0.powi(3) / (zeta3t0 * (1.0 - zeta3t0)).powi(2)
                            * ((6.0 - 4.0 * zeta3t0) / zeta3t0
                                - (5.0 - zeta3t0) / (1.0 - zeta3t0))
                        + 6.0 * zeta2t0.powi(3) / zeta3t0.powi(4) * (1.0 - zeta3t0).ln()))
            + self.sum_xm * (zeta3t2 / (1.0 - zeta3t0) + (zeta3t1 / (1.0 - zeta3t0)).powi(2))
    }
}
/// HsPure
pub struct HsPure {
    sum_xm: f64,
    eta0: (f64, f64),
    eta1: (f64, f64),
    eta2: (f64, f64),
    eta3: (f64, f64),
    eta4: (f64, f64),
}
impl HsPure {
    pub fn mu_k<'a>(
        &mut self,
        eta: f64,
        rho_num: f64,
        m_k: &'a [f64],
        eta_k: &'a [f64],
    ) -> impl Iterator<Item = f64> + use<'a> {
        let coef_m = self.eta0(eta);
        let coef_eta = rho_num * self.sum_xm * self.eta1(eta);
        m_k.iter()
            .zip(eta_k)
            .map(move |(m_k, eta_k)| coef_m * m_k + coef_eta * eta_k)
    }
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
impl HsPure {
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
