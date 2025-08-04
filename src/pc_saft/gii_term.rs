/// GiiTerm
pub struct GiiTerm {
    neg_xm1: Vec<f64>,
    t0d0: (f64, Vec<f64>),
    t0d1: (f64, Vec<f64>),
    t0d2: (f64, Vec<f64>),
    t0d3: (f64, Vec<f64>),
    t1d0: (f64, Vec<f64>),
    t1d1: (f64, Vec<f64>),
    t2d0: (f64, Vec<f64>),
}
impl GiiTerm {
    #[allow(clippy::too_many_arguments)]
    pub fn lngii_mu_k<'a>(
        &mut self,
        zeta2: f64,
        zeta3: f64,
        di: &[f64],
        rho_num: f64,
        m1_k: &'a [f64],
        zeta2_k: &'a [f64],
        zeta3_k: &'a [f64],
    ) -> impl Iterator<Item = f64> + use<'a> {
        self.t0d0_flash(zeta2, zeta3, di);
        let coef_di1 = 1.5 / (1.0 - zeta3).powi(2);
        let coef_di2 = zeta2 / (1.0 - zeta3).powi(3);
        let coef_zeta2 = rho_num
            * self
                .neg_xm1
                .iter()
                .zip(&self.t0d0.1)
                .zip(di)
                .map(|((neg_xm1, t0d0), di)| {
                    neg_xm1 / t0d0 * (di * coef_di1 + di.powi(2) * coef_di2)
                })
                .sum::<f64>();
        let coef_di0 = (1.0 - zeta3).powi(2).recip();
        let coef_di1 = 3.0 * zeta2 / (1.0 - zeta3).powi(3);
        let coef_di2 = 1.5 * zeta2.powi(2) / (1.0 - zeta3).powi(4);
        let coef_zeta3 = rho_num
            * self
                .neg_xm1
                .iter()
                .zip(&self.t0d0.1)
                .zip(di)
                .map(|((neg_xm1, t0d0), di)| {
                    neg_xm1 / t0d0 * (coef_di0 + di * coef_di1 + di.powi(2) * coef_di2)
                })
                .sum::<f64>();
        self.t0d0
            .1
            .clone()
            .into_iter()
            .zip(m1_k)
            .zip(zeta2_k)
            .zip(zeta3_k)
            .map(move |(((t0d0, m1_k), zeta2_k), zeta3_k)| {
                -m1_k * t0d0.ln() + coef_zeta2 * zeta2_k + coef_zeta3 * zeta3_k
            })
    }
    pub fn lngii_t0d0(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) -> f64 {
        self.t0d0_flash(zeta2t0, zeta3t0, dit0);
        self.neg_xm1
            .iter()
            .zip(&self.t0d0.1)
            .map(|(coef, t0d0)| coef * t0d0.ln())
            .sum()
    }
    pub fn lngii_t0d1(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) -> f64 {
        self.t0d0_flash(zeta2t0, zeta3t0, dit0);
        self.t0d1_flash(zeta2t0, zeta3t0, dit0);
        self.neg_xm1
            .iter()
            .zip(&self.t0d0.1)
            .zip(&self.t0d1.1)
            .map(|((coef, t0d0), t0d1)| coef * t0d1 / t0d0)
            .sum()
    }
    pub fn lngii_t0d2(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) -> f64 {
        self.t0d0_flash(zeta2t0, zeta3t0, dit0);
        self.t0d1_flash(zeta2t0, zeta3t0, dit0);
        self.t0d2_flash(zeta2t0, zeta3t0, dit0);
        self.neg_xm1
            .iter()
            .zip(&self.t0d0.1)
            .zip(&self.t0d1.1)
            .zip(&self.t0d2.1)
            .map(|(((coef, t0d0), t0d1), t0d2)| coef * (t0d2 / t0d0 - (t0d1 / t0d0).powi(2)))
            .sum()
    }
    pub fn lngii_t0d3(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) -> f64 {
        self.t0d0_flash(zeta2t0, zeta3t0, dit0);
        self.t0d1_flash(zeta2t0, zeta3t0, dit0);
        self.t0d2_flash(zeta2t0, zeta3t0, dit0);
        self.t0d3_flash(zeta2t0, zeta3t0, dit0);
        self.neg_xm1
            .iter()
            .zip(&self.t0d0.1)
            .zip(&self.t0d1.1)
            .zip(&self.t0d2.1)
            .zip(&self.t0d3.1)
            .map(|((((coef, t0d0), t0d1), t0d2), t0d3)| {
                coef * (t0d3 / t0d0 - 3.0 * t0d2 * t0d1 / t0d0.powi(2)
                    + 2.0 * (t0d1 / t0d0).powi(3))
            })
            .sum()
    }
    pub fn lngii_t1d0(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (&[f64], &[f64]),
    ) -> f64 {
        self.t0d0_flash(zeta2t0, zeta3t0, dit0);
        self.t1d0_flash((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1));
        self.neg_xm1
            .iter()
            .zip(&self.t0d0.1)
            .zip(&self.t1d0.1)
            .map(|((coef, t0d0), t1d0)| coef * t1d0 / t0d0)
            .sum()
    }
    pub fn lngii_t1d1(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (&[f64], &[f64]),
    ) -> f64 {
        self.t0d0_flash(zeta2t0, zeta3t0, dit0);
        self.t0d1_flash(zeta2t0, zeta3t0, dit0);
        self.t1d0_flash((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1));
        self.t1d1_flash((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1));
        self.neg_xm1
            .iter()
            .zip(&self.t0d0.1)
            .zip(&self.t0d1.1)
            .zip(&self.t1d0.1)
            .zip(&self.t1d1.1)
            .map(|((((coef, t0d0), t0d1), t1d0), t1d1)| {
                coef * (t1d1 / t0d0 - t0d1 * t1d0 / t0d0.powi(2))
            })
            .sum()
    }
    pub fn lngii_t2d0(
        &mut self,
        (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
        (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
        (dit0, dit1, dit2): (&[f64], &[f64], &[f64]),
    ) -> f64 {
        self.t0d0_flash(zeta2t0, zeta3t0, dit0);
        self.t1d0_flash((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1));
        self.t2d0_flash(
            (zeta2t0, zeta2t1, zeta2t2),
            (zeta3t0, zeta3t1, zeta3t2),
            (dit0, dit1, dit2),
        );
        self.neg_xm1
            .iter()
            .zip(&self.t0d0.1)
            .zip(&self.t1d0.1)
            .zip(&self.t2d0.1)
            .map(|(((coef, t0d0), t1d0), t2d0)| coef * (t2d0 / t0d0 - (t1d0 / t0d0).powi(2)))
            .sum()
    }
}
impl GiiTerm {
    pub fn new(vec_xm1: &[f64]) -> Self {
        Self {
            neg_xm1: vec_xm1.iter().map(|val| -val).collect(),
            t0d0: (0.0, Vec::new()),
            t0d1: (0.0, Vec::new()),
            t0d2: (0.0, Vec::new()),
            t0d3: (0.0, Vec::new()),
            t1d0: (0.0, Vec::new()),
            t1d1: (0.0, Vec::new()),
            t2d0: (0.0, Vec::new()),
        }
    }
    fn t0d0_flash(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) {
        if zeta3t0 != self.t0d0.0 {
            let coef_di0 = 1.0 / (1.0 - zeta3t0);
            let coef_di1 = 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(2);
            let coef_di2 = 0.5 * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(3);
            self.t0d0 = (
                zeta3t0,
                dit0.iter()
                    .map(|di| coef_di0 + coef_di1 * di + coef_di2 * di.powi(2))
                    .collect(),
            )
        }
    }
    fn t0d1_flash(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) {
        if zeta3t0 != self.t0d1.0 {
            let coef_di0 = zeta3t0 / (1.0 - zeta3t0).powi(2);
            let coef_di1 = 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0);
            let coef_di2 = 0.5 * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(4) * (2.0 + zeta3t0);
            self.t0d1 = (
                zeta3t0,
                dit0.iter()
                    .map(|di| coef_di0 + coef_di1 * di + coef_di2 * di.powi(2))
                    .collect(),
            )
        }
    }
    fn t0d2_flash(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) {
        if zeta3t0 != self.t0d2.0 {
            let coef_di0 = 2.0 * zeta3t0.powi(2) / (1.0 - zeta3t0).powi(3);
            let coef_di1 = 3.0 * zeta2t0 * zeta3t0 / (1.0 - zeta3t0).powi(4) * (2.0 + zeta3t0);
            let coef_di2 =
                zeta2t0.powi(2) / (1.0 - zeta3t0).powi(5) * (1.0 + 4.0 * zeta3t0 + zeta3t0.powi(2));
            self.t0d2 = (
                zeta3t0,
                dit0.iter()
                    .map(|di| coef_di0 + coef_di1 * di + coef_di2 * di.powi(2))
                    .collect(),
            )
        }
    }
    fn t0d3_flash(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: &[f64]) {
        if zeta3t0 != self.t0d3.0 {
            let coef_di0 = 6.0 * zeta3t0.powi(3) / (1.0 - zeta3t0).powi(4);
            let coef_di1 =
                9.0 * zeta2t0 * zeta3t0.powi(2) / (1.0 - zeta3t0).powi(5) * (3.0 + zeta3t0);
            let coef_di2 = 3.0 * zeta2t0.powi(2) * zeta3t0 / (1.0 - zeta3t0).powi(6)
                * (3.0 + 6.0 * zeta3t0 + zeta3t0.powi(2));
            self.t0d3 = (
                zeta3t0,
                dit0.iter()
                    .map(|di| coef_di0 + coef_di1 * di + coef_di2 * di.powi(2))
                    .collect(),
            )
        }
    }
    fn t1d0_flash(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (&[f64], &[f64]),
    ) {
        if zeta3t0 != self.t1d0.0 {
            let coef_di0_dit0 = zeta3t1 / (1.0 - zeta3t0).powi(2);
            let coef_di0_dit1 = 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(2);
            let coef_di1_dit0 = 1.5 * zeta2t1 / (1.0 - zeta3t0).powi(2)
                + 3.0 * zeta2t0 * zeta3t1 / (1.0 - zeta3t0).powi(3);
            let coef_di1_dit1 = zeta2t0.powi(2) / (1.0 - zeta3t0).powi(3);
            let coef_di2_dit0 = zeta2t0 * zeta2t1 / (1.0 - zeta3t0).powi(3)
                + 1.5 * zeta2t0.powi(2) * zeta3t1 / (1.0 - zeta3t0).powi(4);
            self.t1d0 = (
                zeta3t0,
                dit0.iter()
                    .zip(dit1)
                    .map(|(dit0, dit1)| {
                        coef_di0_dit0
                            + coef_di0_dit1 * dit1
                            + coef_di1_dit0 * dit0
                            + coef_di1_dit1 * dit0 * dit1
                            + coef_di2_dit0 * dit0.powi(2)
                    })
                    .collect(),
            )
        }
    }
    fn t1d1_flash(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (&[f64], &[f64]),
    ) {
        if zeta3t0 != self.t1d1.0 {
            let coef_di0_dit0 = zeta3t1 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0);
            let coef_di0_dit1 = 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0);
            let coef_di1_dit0 = 1.5 * zeta2t1 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0)
                + 3.0 * zeta2t0 * zeta3t1 / (1.0 - zeta3t0).powi(4) * (2.0 + zeta3t0);
            let coef_di1_dit1 = zeta2t0.powi(2) / (1.0 - zeta3t0).powi(4) * (2.0 + zeta3t0);
            let coef_di2_dit0 = zeta2t0 * zeta2t1 / (1.0 - zeta3t0).powi(4) * (2.0 + zeta3t0)
                + 1.5 * zeta2t0.powi(2) * zeta3t1 / (1.0 - zeta3t0).powi(5) * (3.0 + zeta3t0);
            self.t1d1 = (
                zeta3t0,
                dit0.iter()
                    .zip(dit1)
                    .map(|(dit0, dit1)| {
                        coef_di0_dit0
                            + coef_di0_dit1 * dit1
                            + coef_di1_dit0 * dit0
                            + coef_di1_dit1 * dit0 * dit1
                            + coef_di2_dit0 * dit0.powi(2)
                    })
                    .collect(),
            )
        }
    }
    fn t2d0_flash(
        &mut self,
        (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
        (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
        (dit0, dit1, dit2): (&[f64], &[f64], &[f64]),
    ) {
        if zeta3t0 != self.t2d0.0 {
            let coef_di0_dit0_ditt0 =
                zeta3t2 / (1.0 - zeta3t0).powi(2) + 2.0 * zeta3t1.powi(2) / (1.0 - zeta3t0).powi(3);
            let coef_di0_dit0_ditt1 = 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(2);
            let coef_di0_dit1_ditt0 = 3.0 * zeta2t1 / (1.0 - zeta3t0).powi(2)
                + 6.0 * zeta2t0 * zeta3t1 / (1.0 - zeta3t0).powi(3);
            let coef_di1_dit0_ditt0 = 1.5 * zeta2t2 / (1.0 - zeta3t0).powi(2)
                + 6.0 * zeta2t1 * zeta3t1 / (1.0 - zeta3t0).powi(3)
                + 3.0 * zeta2t0 * zeta3t2 / (1.0 - zeta3t0).powi(3)
                + 9.0 * zeta2t0 * zeta3t1.powi(2) / (1.0 - zeta3t0).powi(4);
            let coef_di1_dit1_ditt0 = 4.0 * zeta2t0 * zeta2t1 / (1.0 - zeta3t0).powi(3)
                + 6.0 * zeta2t0.powi(2) * zeta3t1 / (1.0 - zeta3t0).powi(4);
            let coef_di0_dit2_ditt0 = zeta2t0.powi(2) / (1.0 - zeta3t0).powi(3);
            let coef_di2_dit0_ditt0 = zeta2t1.powi(2) / (1.0 - zeta3t0).powi(3)
                + zeta2t0 * zeta2t2 / (1.0 - zeta3t0).powi(3)
                + 6.0 * zeta2t0 * zeta2t1 * zeta3t1 / (1.0 - zeta3t0).powi(4)
                + 1.5 * zeta2t0.powi(2) * zeta3t2 / (1.0 - zeta3t0).powi(4)
                + 6.0 * zeta2t0.powi(2) * zeta3t1.powi(2) / (1.0 - zeta3t0).powi(5);
            self.t2d0 = (
                zeta3t0,
                dit0.iter()
                    .zip(dit1)
                    .zip(dit2)
                    .map(|((dit0, dit1), dit2)| {
                        coef_di0_dit0_ditt0
                            + coef_di0_dit0_ditt1 * dit2
                            + coef_di0_dit1_ditt0 * dit1
                            + coef_di1_dit0_ditt0 * dit0
                            + coef_di1_dit1_ditt0 * dit0 * dit1
                            + coef_di0_dit2_ditt0 * (dit1.powi(2) + dit0 * dit2)
                            + coef_di2_dit0_ditt0 * dit0.powi(2)
                    })
                    .collect(),
            )
        }
    }
}
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
