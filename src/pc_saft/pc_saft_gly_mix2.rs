use super::{AssocGlyTerm, AssocType};
use super::{DispTerm, GiiTerm, HsTerm, PolarTerm};
use super::{FRAC_NA_1E30, FRAC_RE30_NA, R};
use super::{PcSaftErr, PcSaftGlyPure};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6};
/// PC-SAFT EOS :: PcSaftGlyMix2
/// ```
/// use thermolib::PcSaftGlyMix2;
/// let mut fluids = PcSaftGlyMix2::new_fluid(
///     [0.5, 0.5],
///     [2.0729, 1.5255],
///     [2.7852, 3.2300],
///     [169.21, 188.90],
///     0.0,
/// ); // CO2+CH3OH_2B
/// fluids.set_2B_assoc_term(1, 0.035176, 2899.5, 1.0, 1.0, 1.0);
/// fluids.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluids.rho().unwrap().round(), 45.0);
/// let py = fluids.tx_flash(298.15, 0.5).unwrap();
/// assert_eq!((py.1 * 1e5).round() / 1e5, 0.99649);
/// let px = fluids.ty_flash(298.15, 0.5).unwrap();
/// assert_eq!((px.1 * 1e5).round() / 1e5, 0.00143);
/// let mut fluids = PcSaftGlyMix2::new_fluid(
///     [0.5, 0.5],
///     [1.5131, 1.5255],
///     [3.1869, 3.2300],
///     [163.33, 188.90],
///     0.0,
/// ); // CO2_QQ+CH3OH_2B
/// fluids.set_2B_assoc_term(1, 0.035176, 2899.5, 1.0, 1.0, 1.0);
/// fluids.set_polar_term([4.4, 0.0], [4, 0]);
/// fluids.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluids.rho().unwrap().round(), 45.0);
/// let py = fluids.tx_flash(298.15, 0.5).unwrap();
/// assert_eq!((py.1 * 1e5).round() / 1e5, 0.99624);
/// let px = fluids.ty_flash(298.15, 0.5).unwrap();
/// assert_eq!((px.1 * 1e5).round() / 1e5, 0.00061);
/// let mut fluids = PcSaftGlyMix2::new_fluid(
///     [0.5, 0.5],
///     [2.7447, 1.5255],
///     [3.2742, 3.2300],
///     [232.99, 188.90],
///     0.0,
/// ); // ACETONE_DD+CH3OH_2B
/// fluids.set_2B_assoc_term(1, 0.035176, 2899.5, 1.0, 1.0, 1.0);
/// fluids.set_polar_term([2.88, 0.0], [2, 0]);
/// fluids.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluids.rho().unwrap().round(), 17116.0);
/// let py = fluids.tx_flash(298.15, 0.5).unwrap();
/// assert_eq!((py.1 * 1e5).round() / 1e5, 0.66555);
/// let px = fluids.ty_flash(298.15, 0.5).unwrap();
/// assert_eq!((px.1 * 1e5).round() / 1e5, 0.06152);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
pub struct PcSaftGlyMix2 {
    x: [f64; 2],
    m: [f64; 2],
    m1: [f64; 2],
    sigma: [f64; 2],
    epsilon: [f64; 2],
    hs: HsTerm,                  // HsTerm
    gii: GiiTerm,                // GiiTerm
    disp: DispTerm,              // DispTerm
    assoc: Option<AssocGlyTerm>, // AssocGlyTerm
    n_assoc: usize,
    polar: Option<PolarTerm>, // PolarTerm
    p: [f64; 2],
    n: [i32; 2],
    m2e1s3_coef: [f64; 3],
    m2e2s3_coef: [f64; 3],
    // state
    temp: f64,
    rho_num: f64,
    dt0: (f64, [f64; 2]),
    dt1: (f64, [f64; 2]),
    dt2: (f64, [f64; 2]),
    zeta_t0_coef: (f64, [[f64; 3]; 2]), // zeta1t0,zeta2t0,zeta3t0
    zeta_t1_coef: (f64, [[f64; 3]; 2]), // zeta1t1,zeta2t1,zeta3t1
    zeta_t2_coef: (f64, [[f64; 3]; 2]), // zeta1t2,zeta2t2,zeta3t2
    eta0_coef: (f64, [f64; 2]),
    zeta1_mu_k: (f64, [f64; 2]),
    zeta2_mu_k: (f64, [f64; 2]),
    zeta3_mu_k: (f64, [f64; 2]),
    is_single_phase: bool,
    // critical point
    omega1: Option<[f64; 2]>,
    temp_c: [f64; 2],
    pres_c: [f64; 2],
}
impl PcSaftGlyMix2 {
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
        let m1 = [m[0] - 1.0, m[1] - 1.0];
        Self {
            x,
            m,
            m1,
            sigma,
            epsilon,
            hs: HsTerm::new(x[0] * m[0] + x[1] * m[1]),
            gii: GiiTerm::new(&[x[0] * m1[0], x[1] * m1[1]]),
            disp: DispTerm::new(
                x[0] * m[0] + x[1] * m[1],
                x[0].powi(2) * m2e1s3_coef[0]
                    + x[1].powi(2) * m2e1s3_coef[1]
                    + x[0] * x[1] * m2e1s3_coef[2],
                x[0].powi(2) * m2e2s3_coef[0]
                    + x[1].powi(2) * m2e2s3_coef[1]
                    + x[0] * x[1] * m2e2s3_coef[2],
            ),
            assoc: None,
            n_assoc: 0,
            polar: None,
            p: [0.0, 0.0],
            n: [0, 0],
            m2e1s3_coef,
            m2e2s3_coef,
            // state
            temp: 0.0,
            rho_num: 0.0,
            dt0: (0.0, [0.0, 0.0]),
            dt1: (0.0, [0.0, 0.0]),
            dt2: (0.0, [0.0, 0.0]),
            zeta_t0_coef: (0.0, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            zeta_t1_coef: (0.0, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            zeta_t2_coef: (0.0, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            eta0_coef: (0.0, [0.0, 0.0]),
            zeta1_mu_k: (0.0, [0.0, 0.0]),
            zeta2_mu_k: (0.0, [0.0, 0.0]),
            zeta3_mu_k: (0.0, [0.0, 0.0]),
            is_single_phase: true,
            // critical point
            omega1: None,
            temp_c: [0.0, 0.0],
            pres_c: [0.0, 0.0],
        }
    }
    fn new_fracs(&self, x: [f64; 2]) -> Self {
        Self {
            x,
            hs: HsTerm::new(x[0] * self.m[0] + x[1] * self.m[1]),
            gii: GiiTerm::new(&[x[0] * self.m1[0], x[1] * self.m1[1]]),
            disp: DispTerm::new(
                x[0] * self.m[0] + x[1] * self.m[1],
                x[0].powi(2) * self.m2e1s3_coef[0]
                    + x[1].powi(2) * self.m2e1s3_coef[1]
                    + x[0] * x[1] * self.m2e1s3_coef[2],
                x[0].powi(2) * self.m2e2s3_coef[0]
                    + x[1].powi(2) * self.m2e2s3_coef[1]
                    + x[0] * x[1] * self.m2e2s3_coef[2],
            ),
            polar: self.polar.as_ref().map(|p| p.new_fracs(&x)),
            assoc: self.assoc.as_ref().map(|a| a.new_x_frac(x[self.n_assoc])),
            ..*self
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)] // For pymethods hhh
impl PcSaftGlyMix2 {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(x: [f64; 2], m: [f64; 2], sigma: [f64; 2], epsilon: [f64; 2], kij: f64) -> Self {
        Self::new_fluid(x, m, sigma, epsilon, kij)
    }
    pub fn set_1_assoc_term(
        &mut self,
        n: usize,
        kappa_AB: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) {
        self.n_assoc = if n == 1 { 1 } else { 0 };
        self.assoc = Some(AssocGlyTerm::new_1_term(
            self.x[self.n_assoc],
            kappa_AB * self.sigma[self.n_assoc].powi(3),
            epsilon_AB,
            c0,
            c1,
            c2,
        ));
    }
    pub fn set_2B_assoc_term(
        &mut self,
        n: usize,
        kappa_AB: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) {
        self.n_assoc = if n == 1 { 1 } else { 0 };
        self.assoc = Some(AssocGlyTerm::new_2B_term(
            self.x[self.n_assoc],
            kappa_AB * self.sigma[self.n_assoc].powi(3),
            epsilon_AB,
            c0,
            c1,
            c2,
        ))
    }
    pub fn set_3B_assoc_term(
        &mut self,
        n: usize,
        kappa_AB: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) {
        self.n_assoc = if n == 1 { 1 } else { 0 };
        self.assoc = Some(AssocGlyTerm::new_3B_term(
            self.x[self.n_assoc],
            kappa_AB * self.sigma[self.n_assoc].powi(3),
            epsilon_AB,
            c0,
            c1,
            c2,
        ))
    }
    pub fn set_4C_assoc_term(
        &mut self,
        n: usize,
        kappa_AB: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) {
        self.n_assoc = if n == 1 { 1 } else { 0 };
        self.assoc = Some(AssocGlyTerm::new_4C_term(
            self.x[self.n_assoc],
            kappa_AB * self.sigma[self.n_assoc].powi(3),
            epsilon_AB,
            c0,
            c1,
            c2,
        ))
    }
    pub fn set_polar_term(&mut self, p: [f64; 2], n: [i32; 2]) {
        self.p = [p[0], p[1]];
        self.n = [n[0], n[1]];
        self.polar = Some(PolarTerm::new(
            &self.x,
            &self.m,
            &self.sigma,
            &self.epsilon,
            &p,
            &n,
        ))
    }
}
fn_tpz_flash_mix2!(PcSaftGlyMix2);
fn_tx_flash_mix2!(PcSaftGlyMix2);
fn_ty_flash_mix2!(PcSaftGlyMix2);
fn_tp_flash!(PcSaftGlyMix2);
fn_single_prop!(PcSaftGlyMix2);
fn_virial_prop!(PcSaftGlyMix2);
impl PcSaftGlyMix2 {
    fn_calc_prop!();
    fn_check_derivatives_mix!();
}
impl PcSaftGlyMix2 {
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
                    [
                        FRAC_PI_6 * self.m[0] * self.dt0.1[0],
                        FRAC_PI_6 * self.m[0] * self.dt0.1[0].powi(2),
                        FRAC_PI_6 * self.m[0] * self.dt0.1[0].powi(3),
                    ],
                    [
                        FRAC_PI_6 * self.m[1] * self.dt0.1[1],
                        FRAC_PI_6 * self.m[1] * self.dt0.1[1].powi(2),
                        FRAC_PI_6 * self.m[1] * self.dt0.1[1].powi(3),
                    ],
                ],
            )
        }
        (
            (self.x[0] * self.zeta_t0_coef.1[0][0] + self.x[1] * self.zeta_t0_coef.1[1][0])
                * rho_num, // zeta1t0
            (self.x[0] * self.zeta_t0_coef.1[0][1] + self.x[1] * self.zeta_t0_coef.1[1][1])
                * rho_num, // zeta2t2
            (self.x[0] * self.zeta_t0_coef.1[0][2] + self.x[1] * self.zeta_t0_coef.1[1][2])
                * rho_num, // zeta3t0
        )
    }
    fn zeta_t1(&mut self, temp: f64, rho_num: f64) -> (f64, f64, f64) {
        if temp != self.zeta_t1_coef.0 {
            self.dt0_flash(temp);
            self.dt1_flash(temp);
            self.zeta_t1_coef = (
                temp,
                [
                    [
                        FRAC_PI_6 * self.m[0] * self.dt1.1[0],
                        FRAC_PI_3 * self.m[0] * self.dt0.1[0] * self.dt1.1[0],
                        FRAC_PI_2 * self.m[0] * self.dt0.1[0].powi(2) * self.dt1.1[0],
                    ],
                    [
                        FRAC_PI_6 * self.m[1] * self.dt1.1[1],
                        FRAC_PI_3 * self.m[1] * self.dt0.1[1] * self.dt1.1[1],
                        FRAC_PI_2 * self.m[1] * self.dt0.1[1].powi(2) * self.dt1.1[1],
                    ],
                ],
            )
        }
        (
            (self.x[0] * self.zeta_t1_coef.1[0][0] + self.x[1] * self.zeta_t1_coef.1[1][0])
                * rho_num, // zeta1t1
            (self.x[0] * self.zeta_t1_coef.1[0][1] + self.x[1] * self.zeta_t1_coef.1[1][1])
                * rho_num, // zeta2t1
            (self.x[0] * self.zeta_t1_coef.1[0][2] + self.x[1] * self.zeta_t1_coef.1[1][2])
                * rho_num, // zeta3t1
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
                    [
                        FRAC_PI_6 * self.m[0] * self.dt2.1[0],
                        FRAC_PI_3
                            * self.m[0]
                            * (self.dt1.1[0].powi(2) + self.dt0.1[0] * self.dt2.1[0]),
                        FRAC_PI_2
                            * self.m[0]
                            * (2.0 * self.dt0.1[0] * self.dt1.1[0].powi(2)
                                + self.dt0.1[0].powi(2) * self.dt2.1[0]),
                    ],
                    [
                        FRAC_PI_6 * self.m[1] * self.dt2.1[1],
                        FRAC_PI_3
                            * self.m[1]
                            * (self.dt1.1[1].powi(2) + self.dt0.1[1] * self.dt2.1[1]),
                        FRAC_PI_2
                            * self.m[1]
                            * (2.0 * self.dt0.1[1] * self.dt1.1[1].powi(2)
                                + self.dt0.1[1].powi(2) * self.dt2.1[1]),
                    ],
                ],
            )
        }
        (
            (self.x[0] * self.zeta_t2_coef.1[0][0] + self.x[1] * self.zeta_t2_coef.1[1][0])
                * rho_num, // zeta1t2
            (self.x[0] * self.zeta_t2_coef.1[0][1] + self.x[1] * self.zeta_t2_coef.1[1][1])
                * rho_num, // zeta2t2
            (self.x[0] * self.zeta_t2_coef.1[0][2] + self.x[1] * self.zeta_t2_coef.1[1][2])
                * rho_num, // zeta3t2
        )
    }
    fn eta0_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta0_coef.0 {
            self.dt0_flash(temp);
            self.eta0_coef = (
                temp,
                [
                    FRAC_PI_6 * self.m[0] * self.dt0.1[0].powi(3),
                    FRAC_PI_6 * self.m[1] * self.dt0.1[1].powi(3),
                ],
            );
        }
        self.x[0] * self.eta0_coef.1[0] + self.x[1] * self.eta0_coef.1[1]
    }
    fn zeta1_mu_k(&mut self, temp: f64) -> [f64; 2] {
        if temp != self.zeta1_mu_k.0 {
            self.dt0_flash(temp);
            self.zeta1_mu_k = (
                temp,
                [
                    FRAC_PI_6 * self.m[0] * self.dt0.1[0],
                    FRAC_PI_6 * self.m[1] * self.dt0.1[1],
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
                    FRAC_PI_6 * self.m[0] * self.dt0.1[0].powi(2),
                    FRAC_PI_6 * self.m[1] * self.dt0.1[1].powi(2),
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
                    FRAC_PI_6 * self.m[0] * self.dt0.1[0].powi(3),
                    FRAC_PI_6 * self.m[1] * self.dt0.1[1].powi(3),
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
                a.t0d0(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[self.n_assoc])
            })
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d0(temp, rho_num, zeta3t0))
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        self.hs.t0d1(rho_num, (zeta1t0, zeta2t0, zeta3t0))
            + self.gii.lngii_t0d1(zeta2t0, zeta3t0, &self.dt0.1)
            + self.disp.t0d1(temp, rho_num, zeta3t0)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t0d1(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[self.n_assoc])
            })
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d1(temp, rho_num, zeta3t0))
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        self.hs.t0d2(rho_num, (zeta1t0, zeta2t0, zeta3t0))
            + self.gii.lngii_t0d2(zeta2t0, zeta3t0, &self.dt0.1)
            + self.disp.t0d2(temp, rho_num, zeta3t0)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t0d2(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[self.n_assoc])
            })
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d2(temp, rho_num, zeta3t0))
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        self.hs.t0d3(rho_num, (zeta1t0, zeta2t0, zeta3t0))
            + self.gii.lngii_t0d3(zeta2t0, zeta3t0, &self.dt0.1)
            + self.disp.t0d3(temp, rho_num, zeta3t0)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t0d3(temp, rho_num, zeta2t0, zeta3t0, self.dt0.1[self.n_assoc])
            })
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d3(temp, rho_num, zeta3t0))
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
                    (self.dt0.1[self.n_assoc], self.dt1.1[self.n_assoc]),
                )
            })
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t1d0(temp, rho_num, zeta3t0, zeta3t1))
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
                    (self.dt0.1[self.n_assoc], self.dt1.1[self.n_assoc]),
                )
            })
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t1d1(temp, rho_num, zeta3t0, zeta3t1))
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
                    (
                        self.dt0.1[self.n_assoc],
                        self.dt1.1[self.n_assoc],
                        self.dt2.1[self.n_assoc],
                    ),
                )
            })
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t2d0(temp, rho_num, zeta3t0, zeta3t1, zeta3t2))
    }
    fn ln_phi(&mut self, temp: f64, pres: f64, eta0_guess: f64) -> [f64; 2] {
        let eta0_coef = self.eta0_coef(temp);
        let rho_num = self.calc_density(temp, pres, eta0_guess / eta0_coef);
        let (zeta1t0, zeta2t0, zeta3t0) = self.zeta_t0(temp, rho_num);
        let zeta1_k = self.zeta1_mu_k(temp);
        let zeta2_k = self.zeta2_mu_k(temp);
        let zeta3_k = self.zeta3_mu_k(temp);
        let m2e1s3_k = [
            2.0 * self.x[0] * self.m2e1s3_coef[0] + self.x[1] * self.m2e1s3_coef[2],
            2.0 * self.x[1] * self.m2e1s3_coef[1] + self.x[0] * self.m2e1s3_coef[2],
        ];
        let m2e2s3_k = [
            2.0 * self.x[0] * self.m2e2s3_coef[0] + self.x[1] * self.m2e2s3_coef[2],
            2.0 * self.x[1] * self.m2e2s3_coef[1] + self.x[0] * self.m2e2s3_coef[2],
        ];
        let t0d1 = self.r_t0d1(temp, rho_num);
        let ln_phi: Vec<f64> = self
            .hs
            .mu_k(
                rho_num,
                (zeta1t0, zeta2t0, zeta3t0),
                &self.m,
                &zeta1_k,
                &zeta2_k,
                &zeta3_k,
            )
            .zip(self.gii.lngii_mu_k(
                zeta2t0,
                zeta3t0,
                &self.dt0.1,
                rho_num,
                &self.m1,
                &zeta2_k,
                &zeta3_k,
            ))
            .zip(self.disp.mu_k(
                temp, rho_num, zeta3t0, &self.m, &zeta3_k, &m2e1s3_k, &m2e2s3_k,
            ))
            .zip(if self.assoc.is_some() {
                self.assoc
                    .as_mut()
                    .unwrap()
                    .mu_k(
                        temp,
                        rho_num,
                        zeta2t0,
                        zeta3t0,
                        self.dt0.1[self.n_assoc],
                        self.n_assoc,
                        &zeta2_k,
                        &zeta3_k,
                    )
                    .collect::<Vec<f64>>()
            } else {
                vec![0.0, 0.0]
            })
            .zip(if self.polar.is_some() {
                self.polar
                    .as_mut()
                    .unwrap()
                    .mu_k(temp, rho_num, zeta3t0, &zeta3_k)
                    .collect::<Vec<f64>>()
            } else {
                vec![0.0, 0.0]
            })
            .map(|((((hs, gii), disp), assoc), polar)| {
                hs + gii + disp + assoc + polar - (1.0 + t0d1).ln()
            })
            .collect();
        [ln_phi[0], ln_phi[1]]
    }
    #[inline]
    fn calc_ln_k(&mut self, temp: f64, pres: f64, x: [f64; 2], y: [f64; 2]) -> [f64; 2] {
        let ln_phi_l = self.new_fracs(x).ln_phi(temp, pres, 0.5);
        let ln_phi_v = self.new_fracs(y).ln_phi(temp, pres, 1e-10);
        [ln_phi_l[0] - ln_phi_v[0], ln_phi_l[1] - ln_phi_v[1]]
    }
    #[allow(non_snake_case)]
    fn guess_ps(&mut self, temp: f64) -> [f64; 2] {
        if self.omega1.is_none() {
            let mut fluid = [
                PcSaftGlyPure::new_fluid(self.m[0], self.sigma[0], self.epsilon[0]),
                PcSaftGlyPure::new_fluid(self.m[1], self.sigma[1], self.epsilon[1]),
            ];
            if let Some(a) = &self.assoc {
                let (assoc_type, kappa_AB_sigma3, epsilon_AB, c0, c1, c2) = a.parameters();
                match assoc_type {
                    AssocType::Type1 => fluid[self.n_assoc].set_1_assoc_term(
                        kappa_AB_sigma3 / self.sigma[self.n_assoc].powi(3),
                        epsilon_AB,
                        c0,
                        c1,
                        c2,
                    ),
                    AssocType::Type2B => fluid[self.n_assoc].set_2B_assoc_term(
                        kappa_AB_sigma3 / self.sigma[self.n_assoc].powi(3),
                        epsilon_AB,
                        c0,
                        c1,
                        c2,
                    ),
                    AssocType::Type3B => fluid[self.n_assoc].set_3B_assoc_term(
                        kappa_AB_sigma3 / self.sigma[self.n_assoc].powi(3),
                        epsilon_AB,
                        c0,
                        c1,
                        c2,
                    ),
                    AssocType::Type4C => fluid[self.n_assoc].set_4C_assoc_term(
                        kappa_AB_sigma3 / self.sigma[self.n_assoc].powi(3),
                        epsilon_AB,
                        c0,
                        c1,
                        c2,
                    ),
                }
            }
            match &self.n[0] {
                2 => fluid[0].set_DD_polar_term(self.p[0]),
                4 => fluid[0].set_QQ_polar_term(self.p[0]),
                _ => (),
            }
            match &self.n[1] {
                2 => fluid[1].set_DD_polar_term(self.p[1]),
                4 => fluid[1].set_QQ_polar_term(self.p[1]),
                _ => (),
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
