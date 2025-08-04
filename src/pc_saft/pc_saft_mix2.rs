use super::{AssocTerm, AssocType};
use super::{DispTerm, GiiTerm, HsTerm};
use super::{FRAC_NA_1E30, FRAC_RE30_NA, R};
use super::{PcSaftErr, PcSaftPure};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::FRAC_PI_6;
#[cfg_attr(feature = "with_pyo3", pyclass)]
pub struct PcSaftMix2 {
    xm: [f64; 2],
    sigma: [f64; 2],
    epsilon: [f64; 2],
    hs: HsTerm,               // HsTerm
    gii: GiiTerm,             // GiiTerm
    disp: DispTerm,           // DispTerm
    assoc: Option<AssocTerm>, // AssocTerm
    m2e1s3_coef: [f64; 3],
    m2e2s3_coef: [f64; 3],
    // state
    temp: f64,
    rho_num: f64,
    dt0: (f64, [f64; 2]),
    dt1: (f64, [f64; 2]),
    dt2: (f64, [f64; 2]),
    zeta_t0_coef: (f64, [f64; 3]), // zeta1t0,zeta2t0,zeta3t0
    zeta_t1_coef: (f64, [f64; 3]), // zeta1t1,zeta2t1,zeta3t1
    zeta_t2_coef: (f64, [f64; 3]), // zeta1t2,zeta2t2,zeta3t2
    eta0_coef: (f64, f64),
    is_single_phase: bool,
    // ln_phi
    m_k: [f64; 2],
    m1_k: [f64; 2],
    m2e1s3_k: [f64; 2],
    m2e2s3_k: [f64; 2],
    zeta1_mu_k: (f64, [f64; 2]),
    zeta2_mu_k: (f64, [f64; 2]),
    zeta3_mu_k: (f64, [f64; 2]),
    // critical point
    omega1: Option<[f64; 2]>,
    temp_c: [f64; 2],
    pres_c: [f64; 2],
}
impl PcSaftMix2 {
    pub fn new_fluid(
        x: [f64; 2],
        m: [f64; 2],
        sigma: [f64; 2],
        epsilon: [f64; 2],
        kij: f64,
    ) -> Self {
        let sigma3 = [sigma[0].powi(3), sigma[1].powi(3)];
        let m2e1s3_coef = [
            m[0].powi(2) * epsilon[0] * sigma3[0],
            m[1].powi(2) * epsilon[1] * sigma3[1],
            2.0 * (m[0] * m[1])
                * ((epsilon[0] * epsilon[1]).sqrt() * (1.0 - kij))
                * ((sigma[0] + sigma[1]).powi(3) / 8.0),
        ];
        let m2e2s3_coef = [
            (m[0] * epsilon[0]).powi(2) * sigma3[0],
            (m[1] * epsilon[1]).powi(2) * sigma3[1],
            2.0 * (m[0] * m[1])
                * (epsilon[0] * epsilon[1] * (1.0 - kij).powi(2))
                * ((sigma[0] + sigma[1]).powi(3) / 8.0),
        ];
        let xm = [x[0] * m[0], x[1] * m[1]];
        let m1 = [m[0] - 1.0, m[1] - 1.0];
        Self {
            xm,
            sigma,
            epsilon,
            hs: HsTerm::new(xm[0] + xm[1]),
            gii: GiiTerm::new(&[x[0] * m1[0], x[1] * m1[1]]),
            disp: DispTerm::new(
                xm[0] + xm[1],
                x[0].powi(2) * m2e1s3_coef[0]
                    + x[1].powi(2) * m2e1s3_coef[1]
                    + x[0] * x[1] * m2e1s3_coef[2],
                x[0].powi(2) * m2e2s3_coef[0]
                    + x[1].powi(2) * m2e2s3_coef[1]
                    + x[0] * x[1] * m2e2s3_coef[2],
            ),
            assoc: None,
            m2e1s3_coef,
            m2e2s3_coef,
            // state
            temp: 0.0,
            rho_num: 0.0,
            dt0: (0.0, [0.0, 0.0]),
            dt1: (0.0, [0.0, 0.0]),
            dt2: (0.0, [0.0, 0.0]),
            zeta_t0_coef: (0.0, [0.0, 0.0, 0.0]),
            zeta_t1_coef: (0.0, [0.0, 0.0, 0.0]),
            zeta_t2_coef: (0.0, [0.0, 0.0, 0.0]),
            eta0_coef: (0.0, 0.0),
            is_single_phase: true,
            // ln_phi
            m_k: m,
            m1_k: m1,
            m2e1s3_k: [
                2.0 * x[0] * m2e1s3_coef[0] + x[1] * m2e1s3_coef[2],
                2.0 * x[1] * m2e1s3_coef[1] + x[0] * m2e1s3_coef[2],
            ],
            m2e2s3_k: [
                2.0 * x[0] * m2e2s3_coef[0] + x[1] * m2e2s3_coef[2],
                2.0 * x[1] * m2e2s3_coef[1] + x[0] * m2e2s3_coef[2],
            ],
            zeta1_mu_k: (0.0, [0.0, 0.0]),
            zeta2_mu_k: (0.0, [0.0, 0.0]),
            zeta3_mu_k: (0.0, [0.0, 0.0]),
            // critical point
            omega1: None,
            temp_c: [0.0, 0.0],
            pres_c: [0.0, 0.0],
        }
    }
    fn new_fracs(&self, x: [f64; 2]) -> Self {
        let xm = [x[0] * self.m_k[0], x[1] * self.m_k[1]];
        Self {
            xm,
            hs: HsTerm::new(xm[0] + xm[1]),
            gii: GiiTerm::new(&[x[0] * self.m1_k[0], x[1] * self.m1_k[1]]),
            disp: DispTerm::new(
                xm[0] + xm[1],
                x[0].powi(2) * self.m2e1s3_coef[0]
                    + x[1].powi(2) * self.m2e1s3_coef[1]
                    + x[0] * x[1] * self.m2e1s3_coef[2],
                x[0].powi(2) * self.m2e2s3_coef[0]
                    + x[1].powi(2) * self.m2e2s3_coef[1]
                    + x[0] * x[1] * self.m2e2s3_coef[2],
            ),
            assoc: self.assoc.as_ref().map(|a| a.new_x_frac(x[0])),
            // state
            zeta_t0_coef: (0.0, [0.0, 0.0, 0.0]),
            zeta_t1_coef: (0.0, [0.0, 0.0, 0.0]),
            zeta_t2_coef: (0.0, [0.0, 0.0, 0.0]),
            eta0_coef: (0.0, 0.0),
            // ln_phi
            m2e1s3_k: [
                2.0 * x[0] * self.m2e1s3_coef[0] + x[1] * self.m2e1s3_coef[2],
                2.0 * x[1] * self.m2e1s3_coef[1] + x[0] * self.m2e1s3_coef[2],
            ],
            m2e2s3_k: [
                2.0 * x[0] * self.m2e2s3_coef[0] + x[1] * self.m2e2s3_coef[2],
                2.0 * x[1] * self.m2e2s3_coef[1] + x[0] * self.m2e2s3_coef[2],
            ],
            ..*self
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)] // For pymethods hhh
impl PcSaftMix2 {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(x: [f64; 2], m: [f64; 2], sigma: [f64; 2], epsilon: [f64; 2], kij: f64) -> Self {
        Self::new_fluid(x, m, sigma, epsilon, kij)
    }
    pub fn set_1_assoc_term(&mut self, x: f64, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_1_term(
            x,
            kappa_AB * self.sigma[0].powi(3),
            epsilon_AB,
        ));
    }
    pub fn set_2B_assoc_term(&mut self, x: f64, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_2B_term(
            x,
            kappa_AB * self.sigma[0].powi(3),
            epsilon_AB,
        ))
    }
    pub fn set_3B_assoc_term(&mut self, x: f64, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_3B_term(
            x,
            kappa_AB * self.sigma[0].powi(3),
            epsilon_AB,
        ))
    }
    pub fn set_4C_assoc_term(&mut self, x: f64, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_4C_term(
            x,
            kappa_AB * self.sigma[0].powi(3),
            epsilon_AB,
        ))
    }
}
fn_tpz_flash_mix2!(PcSaftMix2);
fn_tx_flash_mix2!(PcSaftMix2);
fn_ty_flash_mix2!(PcSaftMix2);
fn_tp_flash!(PcSaftMix2);
fn_single_prop!(PcSaftMix2);
fn_virial_prop!(PcSaftMix2);
impl PcSaftMix2 {
    fn_calc_prop!();
    fn_check_derivatives_mix!();
}
impl PcSaftMix2 {
    fn dt0_flash(&mut self, temp: f64) {
        if temp != self.dt0.0 {
            self.dt0 = (
                temp,
                [
                    self.sigma[0] * (1.0 - 0.12 * (-3.0 * self.epsilon[0] / temp).exp()),
                    self.sigma[1] * (1.0 - 0.12 * (-3.0 * self.epsilon[1] / temp).exp()),
                ],
            )
        }
    }
    fn dt1_flash(&mut self, temp: f64) {
        if temp != self.dt1.0 {
            let epsilon_temp_plus = [3.0 * self.epsilon[0] / temp, 3.0 * self.epsilon[1] / temp];
            self.dt1 = (
                temp,
                [
                    self.sigma[0]
                        * (-0.12 * (-epsilon_temp_plus[0]).exp())
                        * (epsilon_temp_plus[0]),
                    self.sigma[1]
                        * (-0.12 * (-epsilon_temp_plus[1]).exp())
                        * (epsilon_temp_plus[1]),
                ],
            )
        }
    }
    fn dt2_flash(&mut self, temp: f64) {
        if temp != self.dt2.0 {
            let epsilon_temp_plus = [3.0 * self.epsilon[0] / temp, 3.0 * self.epsilon[1] / temp];
            self.dt2 = (
                temp,
                [
                    self.sigma[0]
                        * (-0.12 * (-epsilon_temp_plus[0]).exp())
                        * epsilon_temp_plus[0]
                        * (epsilon_temp_plus[0] - 2.0),
                    self.sigma[1]
                        * (-0.12 * (-epsilon_temp_plus[1]).exp())
                        * epsilon_temp_plus[1]
                        * (epsilon_temp_plus[1] - 2.0),
                ],
            )
        }
    }
    fn zeta_t0(&mut self, temp: f64, rho_num: f64) -> (f64, f64, f64) {
        if temp != self.zeta_t0_coef.0 {
            self.dt0_flash(temp);
            self.zeta_t0_coef = (
                temp,
                [
                    FRAC_PI_6 * (self.xm[0] * self.dt0.1[0] + self.xm[1] * self.dt0.1[1]),
                    FRAC_PI_6
                        * (self.xm[0] * self.dt0.1[0].powi(2) + self.xm[1] * self.dt0.1[1].powi(2)),
                    FRAC_PI_6
                        * (self.xm[0] * self.dt0.1[0].powi(3) + self.xm[1] * self.dt0.1[1].powi(3)),
                ],
            )
        }
        (
            self.zeta_t0_coef.1[0] * rho_num, // zeta1t0
            self.zeta_t0_coef.1[1] * rho_num, // zeta2t0
            self.zeta_t0_coef.1[2] * rho_num, // zeta3t0
        )
    }
    fn zeta_t1(&mut self, temp: f64, rho_num: f64) -> (f64, f64, f64) {
        if temp != self.zeta_t1_coef.0 {
            self.dt0_flash(temp);
            self.dt1_flash(temp);
            self.zeta_t1_coef = (
                temp,
                [
                    FRAC_PI_6 * (self.xm[0] * self.dt1.1[0] + self.xm[1] * self.dt1.1[1]),
                    FRAC_PI_6
                        * (self.xm[0] * self.dt0.1[0] * self.dt1.1[0]
                            + self.xm[1] * self.dt0.1[1] * self.dt1.1[1]),
                    FRAC_PI_6
                        * (self.xm[0] * self.dt0.1[0].powi(2) * self.dt1.1[0]
                            + self.xm[1] * self.dt0.1[1].powi(2) * self.dt1.1[1]),
                ],
            )
        }
        (
            self.zeta_t1_coef.1[0] * rho_num, // zeta1t1
            self.zeta_t1_coef.1[1] * rho_num, // zeta2t1
            self.zeta_t1_coef.1[2] * rho_num, // zeta3t1
        )
    }
    fn zeta_t2(&mut self, temp: f64, rho_num: f64) -> (f64, f64, f64) {
        if temp != self.zeta_t2_coef.0 {
            self.dt0_flash(temp);
            self.dt1_flash(temp);
            self.dt2_flash(temp);
            self.zeta_t2_coef = (
                temp,
                [
                    FRAC_PI_6 * (self.xm[0] * self.dt2.1[0] + self.xm[1] * self.dt2.1[1]),
                    FRAC_PI_6
                        * (self.xm[0] * (self.dt1.1[0].powi(2) + self.dt0.1[0] * self.dt2.1[0])
                            + self.xm[1] * (self.dt1.1[1].powi(2) + self.dt0.1[1] * self.dt2.1[1])),
                    FRAC_PI_6
                        * (self.xm[0]
                            * (2.0 * self.dt0.1[0] * self.dt1.1[0].powi(2)
                                + self.dt0.1[0].powi(2) * self.dt2.1[0])
                            + self.xm[1]
                                * (2.0 * self.dt0.1[1] * self.dt1.1[1].powi(2)
                                    + self.dt0.1[1].powi(2) * self.dt2.1[1])),
                ],
            )
        }
        (
            self.zeta_t2_coef.1[0] * rho_num, // zeta1t2
            self.zeta_t2_coef.1[1] * rho_num, // zeta2t2
            self.zeta_t2_coef.1[2] * rho_num, // zeta3t2
        )
    }
    fn eta0_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta0_coef.0 {
            self.dt0_flash(temp);
            self.eta0_coef = (
                temp,
                FRAC_PI_6
                    * (self.xm[0] * self.dt0.1[0].powi(3) + self.xm[1] * self.dt0.1[1].powi(3)),
            );
        }
        self.eta0_coef.1
    }
    fn zeta1_mu_k(&mut self, temp: f64) -> [f64; 2] {
        if temp != self.zeta1_mu_k.0 {
            self.dt0_flash(temp);
            self.zeta1_mu_k = (
                temp,
                [
                    FRAC_PI_6 * self.m_k[0] * self.dt0.1[0],
                    FRAC_PI_6 * self.m_k[1] * self.dt0.1[1],
                ],
            )
        }
        self.zeta1_mu_k.1
    }
    fn zeta2_mu_k(&mut self, temp: f64) -> [f64; 2] {
        if temp != self.zeta2_mu_k.0 {
            self.dt0_flash(temp);
            self.zeta2_mu_k = (
                temp,
                [
                    FRAC_PI_6 * self.m_k[0] * self.dt0.1[0].powi(2),
                    FRAC_PI_6 * self.m_k[1] * self.dt0.1[1].powi(2),
                ],
            )
        }
        self.zeta2_mu_k.1
    }
    fn zeta3_mu_k(&mut self, temp: f64) -> [f64; 2] {
        if temp != self.zeta3_mu_k.0 {
            self.dt0_flash(temp);
            self.zeta3_mu_k = (
                temp,
                [
                    FRAC_PI_6 * self.m_k[0] * self.dt0.1[0].powi(3),
                    FRAC_PI_6 * self.m_k[1] * self.dt0.1[1].powi(3),
                ],
            )
        }
        self.zeta3_mu_k.1
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        self.hs.t0d0(rho_num, (zeta1t0, zeta2t0, zeta3t0))
            + self.gii.lngii_t0d0(zeta2t0, zeta3t0, &self.dt0.1)
            + self.disp.t0d0(temp, rho_num, zeta3t0)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t0d0(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[0])
            })
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        self.hs.t0d1(rho_num, (zeta1t0, zeta2t0, zeta3t0))
            + self.gii.lngii_t0d1(zeta2t0, zeta3t0, &self.dt0.1)
            + self.disp.t0d1(temp, rho_num, zeta3t0)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t0d1(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[0])
            })
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        self.hs.t0d2(rho_num, (zeta1t0, zeta2t0, zeta3t0))
            + self.gii.lngii_t0d2(zeta2t0, zeta3t0, &self.dt0.1)
            + self.disp.t0d2(temp, rho_num, zeta3t0)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t0d2(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[0])
            })
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        self.hs.t0d3(rho_num, (zeta1t0, zeta2t0, zeta3t0))
            + self.gii.lngii_t0d3(zeta2t0, zeta3t0, &self.dt0.1)
            + self.disp.t0d3(temp, rho_num, zeta3t0)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t0d3(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[0])
            })
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        let (zeta1t1, zeta2t1, zeta3t1) = self.zeta_t1(temp, rho_num);
        self.hs.t1d0(
            rho_num,
            (zeta1t0, zeta2t0, zeta3t0),
            (zeta1t1, zeta2t1, zeta3t1),
        ) + self.gii.lngii_t1d0(
            (zeta2t0, zeta2t1),
            (zeta3t0, zeta3t1),
            (&self.dt0.1, &self.dt1.1),
        ) + self.disp.t1d0(temp, rho_num, zeta3t0, zeta3t1)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t1d0(
                    temp,
                    rho_num,
                    (zeta2t0, zeta2t1),
                    (zeta3t0, zeta3t1),
                    (self.dt0.1[0], self.dt1.1[0]),
                )
            })
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        let (zeta1t1, zeta2t1, zeta3t1) = self.zeta_t1(temp, rho_num);
        self.hs.t1d1(
            rho_num,
            (zeta1t0, zeta2t0, zeta3t0),
            (zeta1t1, zeta2t1, zeta3t1),
        ) + self.gii.lngii_t1d1(
            (zeta2t0, zeta2t1),
            (zeta3t0, zeta3t1),
            (&self.dt0.1, &self.dt1.1),
        ) + self.disp.t1d1(temp, rho_num, zeta3t0, zeta3t1)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t1d1(
                    temp,
                    rho_num,
                    (zeta2t0, zeta2t1),
                    (zeta3t0, zeta3t1),
                    (self.dt0.1[0], self.dt1.1[0]),
                )
            })
    }
    fn r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        let (zeta1t1, zeta2t1, zeta3t1) = self.zeta_t1(temp, rho_num);
        let (zeta1t2, zeta2t2, zeta3t2) = self.zeta_t2(temp, rho_num);
        self.hs.t2d0(
            rho_num,
            (zeta1t0, zeta2t0, zeta3t0),
            (zeta1t1, zeta2t1, zeta3t1),
            (zeta1t2, zeta2t2, zeta3t2),
        ) + self.gii.lngii_t2d0(
            (zeta2t0, zeta2t1, zeta2t2),
            (zeta3t0, zeta3t1, zeta3t2),
            (&self.dt0.1, &self.dt1.1, &self.dt2.1),
        ) + self.disp.t2d0(temp, rho_num, zeta3t0, zeta3t1, zeta3t2)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t2d0(
                    temp,
                    rho_num,
                    (zeta2t0, zeta2t1, zeta2t2),
                    (zeta3t0, zeta3t1, zeta3t2),
                    (self.dt0.1[0], self.dt1.1[0], self.dt2.1[0]),
                )
            })
    }
    fn ln_phi(&mut self, temp: f64, pres: f64, eta0_guess: f64) -> [f64; 2] {
        let eta0_coef = self.eta0_coef(temp);
        let rho_num = self.calc_density(temp, pres, eta0_guess / eta0_coef);
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        let zeta1_k = self.zeta1_mu_k(temp);
        let zeta2_k = self.zeta2_mu_k(temp);
        let zeta3_k = self.zeta3_mu_k(temp);
        let t0d1 = self.r_t0d1(temp, rho_num);
        let ln_phi: Vec<f64> = self
            .hs
            .mu_k(
                rho_num,
                (zeta1t0, zeta2t0, zeta3t0),
                &self.m_k,
                &zeta1_k,
                &zeta2_k,
                &zeta3_k,
            )
            .zip(self.gii.lngii_mu_k(
                zeta2t0,
                zeta3t0,
                &self.dt0.1,
                rho_num,
                &self.m1_k,
                &zeta2_k,
                &zeta3_k,
            ))
            .zip(self.disp.mu_k(
                temp,
                rho_num,
                zeta3t0,
                &self.m_k,
                &zeta3_k,
                &self.m2e1s3_k,
                &self.m2e2s3_k,
            ))
            .map(|((hs, gii), disp)| hs + gii + disp - (1.0 + t0d1).ln())
            .collect();
        [
            ln_phi[0]
                + self.assoc.as_mut().map_or(0.0, |a| {
                    a.mu_k(
                        temp,
                        rho_num,
                        zeta2t0,
                        zeta3t0,
                        self.dt0.1[0],
                        zeta2_k[0],
                        zeta3_k[0],
                    )
                }),
            ln_phi[1],
        ]
    }
    fn calc_ln_k(&mut self, temp: f64, pres: f64, x: [f64; 2], y: [f64; 2]) -> [f64; 2] {
        let ln_phi_l = self.new_fracs(x).ln_phi(temp, pres, 0.5);
        let ln_phi_v = self.new_fracs(y).ln_phi(temp, pres, 1e-10);
        [ln_phi_l[0] - ln_phi_v[0], ln_phi_l[1] - ln_phi_v[1]]
    }
    #[allow(non_snake_case)]
    fn guess_ps(&mut self, temp: f64) -> [f64; 2] {
        if self.omega1.is_none() {
            let mut fluid = [
                PcSaftPure::new_fluid(self.m_k[0], self.sigma[0], self.epsilon[0]),
                PcSaftPure::new_fluid(self.m_k[1], self.sigma[1], self.epsilon[1]),
            ];
            if let Some(a) = &self.assoc {
                let (assoc_type, kappa_AB_sigma3, epsilon_AB) = a.parameters();
                match assoc_type {
                    AssocType::Type1 => fluid[0]
                        .set_1_assoc_term(kappa_AB_sigma3 / self.sigma[0].powi(3), epsilon_AB),
                    AssocType::Type2B => fluid[0]
                        .set_2B_assoc_term(kappa_AB_sigma3 / self.sigma[0].powi(3), epsilon_AB),
                    AssocType::Type3B => fluid[0]
                        .set_3B_assoc_term(kappa_AB_sigma3 / self.sigma[0].powi(3), epsilon_AB),
                    AssocType::Type4C => fluid[0]
                        .set_4C_assoc_term(kappa_AB_sigma3 / self.sigma[0].powi(3), epsilon_AB),
                }
            }
            fluid[0].c_flash().unwrap();
            fluid[1].c_flash().unwrap();
            self.temp_c = [fluid[0].T().unwrap(), fluid[1].T().unwrap()];
            self.pres_c = [fluid[0].p().unwrap(), fluid[1].p().unwrap()];
            fluid[0].t_flash(0.7 * self.temp_c[0]).unwrap();
            fluid[1].t_flash(0.7 * self.temp_c[1]).unwrap();
            self.omega1 = Some([
                -(fluid[0].p_s().unwrap() / self.pres_c[0]).log10(),
                -(fluid[1].p_s().unwrap() / self.pres_c[1]).log10(),
            ]);
        }
        [
            self.pres_c[0]
                * 10_f64.powf(7.0 / 3.0 * self.omega1.unwrap()[0] * (1.0 - self.temp_c[0] / temp)),
            self.pres_c[1]
                * 10_f64.powf(7.0 / 3.0 * self.omega1.unwrap()[1] * (1.0 - self.temp_c[1] / temp)),
        ]
    }
}
#[cfg(test)]
mod tests {
    #[test]
    fn test_pc_saft_mix2() {
        let mut fluids = super::PcSaftMix2::new_fluid(
            [0.5, 0.5],
            [1.0000, 1.6069],
            [3.7039, 3.5206],
            [150.03, 191.42],
            0.0,
        );
        let temp = 199.92;
        // test tpz_flash
        let p_data = [
            3.62e5, 4.42e5, 6.8e5, 10.9e5, 17e5, 23.8e5, 34e5, 40.8e5, 47.65e5, 48.9e5, 49.4e5,
            49.8e5, 50.35e5,
        ];
        [
            (0.0291, 0.3962),
            (0.0452, 0.5041),
            (0.0932, 0.6756),
            (0.1765, 0.7959),
            (0.3015, 0.8682),
            (0.4414, 0.9058),
            (0.6470, 0.9357),
            (0.7751, 0.9488),
            (0.8907, 0.9597),
            (0.9101, 0.9613),
            (0.9178, 0.9619),
            (0.9239, 0.9623),
            (0.9323, 0.9627),
        ]
        .iter()
        .zip(p_data.iter())
        .map(|(&(x, y), &p)| {
            let xy = fluids.tpz_flash(temp, p).unwrap();
            assert_eq!((xy[0] * 1e4).round() / 1e4, x);
            assert_eq!((xy[1] * 1e4).round() / 1e4, y);
        })
        .count();
        // test tx_flash
        let x_data = [
            0.0214, 0.0512, 0.1039, 0.1875, 0.31, 0.4526, 0.6601, 0.7852, 0.8942, 0.9126, 0.9175,
            0.9222, 0.9319,
        ];
        [
            (3.24e5, 0.3258),
            (4.72e5, 0.5351),
            (7.33e5, 0.6986),
            (11.44e5, 0.8054),
            (17.41e5, 0.8712),
            (24.35e5, 0.9079),
            (34.67e5, 0.9371),
            (41.37e5, 0.9498),
            (47.87e5, 0.9600),
            (49.06e5, 0.9615),
            (49.38e5, 0.9619),
            (49.69e5, 0.9622),
            (50.33e5, 0.9627),
        ]
        .iter()
        .zip(x_data.iter())
        .map(|(&(p, y), &x)| {
            let py = fluids.tx_flash(temp, x).unwrap();
            assert_eq!((py[0] / 1e3).round() * 1e3, p);
            assert_eq!((py[1] * 1e4).round() / 1e4, y);
        })
        .count();
        // test ty_flash
        let y_data = [
            0.3005, 0.5098, 0.68, 0.7957, 0.8679, 0.9052, 0.9337, 0.9461, 0.9562, 0.9584, 0.9578,
            0.9575, 0.9577,
        ];
        [
            (3.12e5, 0.0190),
            (4.47e5, 0.0462),
            (6.90e5, 0.0951),
            (10.89e5, 0.1763),
            (16.97e5, 0.3008),
            (23.65e5, 0.4384),
            (33.12e5, 0.6296),
            (39.28e5, 0.7474),
        ]
        .iter()
        .zip(y_data.iter())
        .map(|(&(p, x), &y)| {
            let px = fluids.ty_flash(temp, y).unwrap();
            assert_eq!((px[0] / 1e3).round() * 1e3, p);
            assert_eq!((px[1] * 1e4).round() / 1e4, x);
        })
        .count();
    }
}
