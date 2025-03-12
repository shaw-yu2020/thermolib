use super::PcSaftErr;
use super::{DispTerm, GiiTerm, HsTerm};
use super::{FRAC_NA_1E30, FRAC_RE30_NA, R};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::FRAC_PI_6;
/// PC-SAFT EOS
/// ```
/// use thermolib::PcSaftGlyPure;
/// let mut methanol = PcSaftGlyPure::new_fluid(1.5255, 3.23, 188.9);
/// methanol.set_2B_assoc_term(0.035176, 2899.5); // kappa_AB epsilon_AB
/// methanol.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(methanol.rho().unwrap().round(), 24676.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)] // For pyclass hhh
pub struct PcSaftGlyPure {
    m: f64,
    sigma3: f64,
    epsilon: f64,
    hs: HsTerm,     // HsTerm
    gii: GiiTerm,   // GiiTerm
    disp: DispTerm, // DispTerm
    eta0: f64,
    eta1: f64,
    eta2: f64,
    // state
    temp: f64,
    rho_num: f64,
    rhov_num: f64,
    rhol_num: f64,
    is_single_phase: bool,
    // modified aly_lee_cv0
    // for fn_aly_lee_cp0!
    cv_B: f64, // cv_B=B/R-1
    cv_C: f64, // cv_C=C/R
    cv_D: f64, // cv_D=D
    cv_E: f64, // cv_E=E/R
    cv_F: f64, // cv_F=F
    // parameters for gly
    c0: f64,
    c1: f64,
    c2: f64,
    // AssocGlyTerm
    assoc_type: AssocType,
    kappa_AB_sigma3: f64,
    epsilon_AB_temp: f64,
    epsilon_AB: f64,
    XA: f64,
    // cached variables
    xt1: (f64, f64),
    xt2: (f64, f64),
    xt3: (f64, f64),
    xt4: (f64, f64),
    x_t0d1: (f64, f64),
    x_t0d2: (f64, f64),
    x_t0d3: (f64, f64),
    x_t0d4: (f64, f64),
    x_t1d0: (f64, f64),
    x_t1d1: (f64, f64),
    x_t1d2: (f64, f64),
    x_t1d3: (f64, f64),
    x_t2d0: (f64, f64),
    t_t0d0: (f64, f64),
    t_t0d1: (f64, f64),
    t_t0d2: (f64, f64),
    t_t0d3: (f64, f64),
    t_t0d4: (f64, f64),
    t_t1d0: (f64, f64),
    t_t1d1: (f64, f64),
    t_t1d2: (f64, f64),
    t_t1d3: (f64, f64),
    t_t2d0: (f64, f64),
    g_t0d0: (f64, f64),
    g_t0d1: (f64, f64),
    g_t0d2: (f64, f64),
    g_t0d3: (f64, f64),
    g_t0d4: (f64, f64),
    g_t1d0: (f64, f64),
    g_t1d1: (f64, f64),
    g_t1d2: (f64, f64),
    g_t1d3: (f64, f64),
    g_t2d0: (f64, f64),
}
impl PcSaftGlyPure {
    pub fn new_fluid(m: f64, sigma: f64, epsilon: f64) -> Self {
        let sigma3 = sigma.powi(3);
        Self {
            m,
            sigma3,
            epsilon,
            hs: HsTerm::new(),
            gii: GiiTerm::new(),
            disp: DispTerm::new(m, sigma3, epsilon),
            eta0: 0.0,
            eta1: 0.0,
            eta2: 0.0,
            // state
            temp: 1.0,
            rho_num: 1E-10,
            rhov_num: 0.0,
            rhol_num: 0.0,
            is_single_phase: true,
            // modified aly_lee_cv0
            // for fn_aly_lee_cp0!
            cv_B: 1.5,
            cv_C: 0.0,
            cv_D: 1.0,
            cv_E: 0.0,
            cv_F: 1.0,
            // parameters for gly
            c0: 1.0,
            c1: -0.5,
            c2: 0.0,
            // AssocGlyTerm
            assoc_type: AssocType::Type0,
            kappa_AB_sigma3: 0.0,
            epsilon_AB_temp: 0.0,
            epsilon_AB: 0.0,
            XA: 1.0,
            // cached variables
            xt1: (0.0, 0.0),
            xt2: (0.0, 0.0),
            xt3: (0.0, 0.0),
            xt4: (0.0, 0.0),
            x_t0d1: (0.0, 0.0),
            x_t0d2: (0.0, 0.0),
            x_t0d3: (0.0, 0.0),
            x_t0d4: (0.0, 0.0),
            x_t1d0: (0.0, 0.0),
            x_t1d1: (0.0, 0.0),
            x_t1d2: (0.0, 0.0),
            x_t1d3: (0.0, 0.0),
            x_t2d0: (0.0, 0.0),
            t_t0d0: (0.0, 0.0),
            t_t0d1: (0.0, 0.0),
            t_t0d2: (0.0, 0.0),
            t_t0d3: (0.0, 0.0),
            t_t0d4: (0.0, 0.0),
            t_t1d0: (0.0, 0.0),
            t_t1d1: (0.0, 0.0),
            t_t1d2: (0.0, 0.0),
            t_t1d3: (0.0, 0.0),
            t_t2d0: (0.0, 0.0),
            g_t0d0: (0.0, 0.0),
            g_t0d1: (0.0, 0.0),
            g_t0d2: (0.0, 0.0),
            g_t0d3: (0.0, 0.0),
            g_t0d4: (0.0, 0.0),
            g_t1d0: (0.0, 0.0),
            g_t1d1: (0.0, 0.0),
            g_t1d2: (0.0, 0.0),
            g_t1d3: (0.0, 0.0),
            g_t2d0: (0.0, 0.0),
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)] // For pymethods hhh
impl PcSaftGlyPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_2B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type2B;
        self.kappa_AB_sigma3 = kappa_AB * self.sigma3;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_3B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type3B;
        self.kappa_AB_sigma3 = kappa_AB * self.sigma3;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_gly_params(&mut self, c0: f64, c1: f64, c2: f64) {
        self.c0 = c0;
        self.c1 = -2.0 * c0 + 1.5 * c1;
        self.c2 = c0 - 1.5 * c1 + 0.5 * c2;
    }
}
fn_vec!(PcSaftGlyPure);
fn_flash!(PcSaftGlyPure);
fn_c_flash!(PcSaftGlyPure);
fn_single_prop!(PcSaftGlyPure);
fn_double_prop!(PcSaftGlyPure);
fn_virial_prop!(PcSaftGlyPure);
fn_aly_lee_cp0!(PcSaftGlyPure);
impl PcSaftGlyPure {
    fn_calc_prop!();
    fn_check_derivatives!();
}
impl PcSaftGlyPure {
    fn eta_flash(&mut self, temp: f64, rho_num: f64) {
        if temp != self.temp || rho_num != self.rho_num {
            (self.temp, self.rho_num) = (temp, rho_num);
            let epsilon_temp = self.epsilon / temp;
            let d = 1.0 - 0.12 * (-3.0 * epsilon_temp).exp();
            let d1 = -0.36 * (-3.0 * epsilon_temp).exp() * epsilon_temp;
            let d2 = d1 * (3.0 * epsilon_temp - 2.0);
            let rho_plus = FRAC_PI_6 * rho_num * self.m * self.sigma3;
            self.eta0 = rho_plus * d.powi(3);
            self.eta1 = rho_plus * 3.0 * d.powi(2) * d1;
            self.eta2 = rho_plus * (6.0 * d * d1.powi(2) + 3.0 * d.powi(2) * d2);
            match self.assoc_type {
                AssocType::Type0 => (),
                AssocType::Type2B => {
                    self.epsilon_AB_temp = self.epsilon_AB / temp;
                    self.XA = (-1.0 + (1.0 + 4.0 * self.t_t0d0()).sqrt()) / (2.0 * self.t_t0d0());
                }
                AssocType::Type3B => {
                    self.epsilon_AB_temp = self.epsilon_AB / temp;
                    self.XA = (-(1.0 - self.t_t0d0())
                        + ((1.0 - self.t_t0d0()).powi(2) + 8.0 * self.t_t0d0()).sqrt())
                        / (4.0 * self.t_t0d0());
                }
            }
        }
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d0(self.eta0)
            + (1.0 - self.m) * self.gii.lngii_t0d0(self.eta0)
            + self.disp.t0d0(temp, rho_num, self.eta0)
            + self.assoc_t0d0()
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d1(self.eta0)
            + (1.0 - self.m) * self.gii.lngii_t0d1(self.eta0)
            + self.disp.t0d1(temp, rho_num, self.eta0)
            + self.assoc_t0d1()
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d2(self.eta0)
            + (1.0 - self.m) * self.gii.lngii_t0d2(self.eta0)
            + self.disp.t0d2(temp, rho_num, self.eta0)
            + self.assoc_t0d2()
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d3(self.eta0)
            + (1.0 - self.m) * self.gii.lngii_t0d3(self.eta0)
            + self.disp.t0d3(temp, rho_num, self.eta0)
            + self.assoc_t0d3()
    }
    fn r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d4(self.eta0)
            + (1.0 - self.m) * self.gii.lngii_t0d4(self.eta0)
            + self.disp.t0d4(temp, rho_num, self.eta0)
            + self.assoc_t0d4()
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d0(self.eta0, self.eta1)
            + (1.0 - self.m) * self.gii.lngii_t1d0(self.eta0, self.eta1)
            + self.disp.t1d0(temp, rho_num, self.eta0, self.eta1)
            + self.assoc_t1d0()
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d1(self.eta0, self.eta1)
            + (1.0 - self.m) * self.gii.lngii_t1d1(self.eta0, self.eta1)
            + self.disp.t1d1(temp, rho_num, self.eta0, self.eta1)
            + self.assoc_t1d1()
    }
    fn r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d2(self.eta0, self.eta1)
            + (1.0 - self.m) * self.gii.lngii_t1d2(self.eta0, self.eta1)
            + self.disp.t1d2(temp, rho_num, self.eta0, self.eta1)
            + self.assoc_t1d2()
    }
    fn r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d3(self.eta0, self.eta1)
            + (1.0 - self.m) * self.gii.lngii_t1d3(self.eta0, self.eta1)
            + self.disp.t1d3(temp, rho_num, self.eta0, self.eta1)
            + self.assoc_t1d3()
    }
    fn r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t2d0(self.eta0, self.eta1, self.eta2)
            + (1.0 - self.m) * self.gii.lngii_t2d0(self.eta0, self.eta1, self.eta2)
            + self
                .disp
                .t2d0(temp, rho_num, self.eta0, self.eta1, self.eta2)
            + self.assoc_t2d0()
    }
}
enum AssocType {
    Type0,
    Type2B,
    Type3B,
}
impl PcSaftGlyPure {
    fn assoc_t0d0(&self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.XA.ln() - self.XA + 1.0,
            AssocType::Type3B => {
                2.0 * self.XA.ln() + (2.0 * self.XA - 1.0).ln() - 2.0 * self.XA + 2.0
            }
        }
    }
    fn assoc_t0d1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d1::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d1::<1>(self.XA) + self.site_t0d1::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t0d2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d2::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d2::<1>(self.XA) + self.site_t0d2::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t0d3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d3::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d3::<1>(self.XA) + self.site_t0d3::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t0d4(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d4::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d4::<1>(self.XA) + self.site_t0d4::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d0::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d0::<1>(self.XA) + self.site_t1d0::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d1::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d1::<1>(self.XA) + self.site_t1d1::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d2::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d2::<1>(self.XA) + self.site_t1d2::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d3::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d3::<1>(self.XA) + self.site_t1d3::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t2d0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t2d0::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t2d0::<1>(self.XA) + self.site_t2d0::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
}
impl PcSaftGlyPure {
    fn site_t0d1<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x1::<C>(x) * self.x_t0d1()
    }
    fn site_t0d2<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x2::<C>(x) * self.x_t0d1().powi(2) + self.site_x1::<C>(x) * self.x_t0d2()
    }
    fn site_t0d3<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x3::<C>(x) * self.x_t0d1().powi(3)
            + 3.0 * self.site_x2::<C>(x) * self.x_t0d1() * self.x_t0d2()
            + self.site_x1::<C>(x) * self.x_t0d3()
    }
    fn site_t0d4<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x4::<C>(x) * self.x_t0d1().powi(4)
            + 6.0 * self.site_x3::<C>(x) * self.x_t0d1().powi(2) * self.x_t0d2()
            + 3.0 * self.site_x2::<C>(x) * self.x_t0d2().powi(2)
            + 4.0 * self.site_x2::<C>(x) * self.x_t0d1() * self.x_t0d3()
            + self.site_x1::<C>(x) * self.x_t0d4()
    }
    fn site_t1d0<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x1::<C>(x) * self.x_t1d0()
    }
    fn site_t1d1<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x2::<C>(x) * self.x_t1d0() * self.x_t0d1() + self.site_x1::<C>(x) * self.x_t1d1()
    }
    fn site_t1d2<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x3::<C>(x) * self.x_t1d0() * self.x_t0d1().powi(2)
            + 2.0 * self.site_x2::<C>(x) * self.x_t1d1() * self.x_t0d1()
            + self.site_x2::<C>(x) * self.x_t1d0() * self.x_t0d2()
            + self.site_x1::<C>(x) * self.x_t1d2()
    }
    fn site_t1d3<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x4::<C>(x) * self.x_t1d0() * self.x_t0d1().powi(3)
            + 3.0 * self.site_x3::<C>(x) * self.x_t1d1() * self.x_t0d1().powi(2)
            + 3.0 * self.site_x3::<C>(x) * self.x_t1d0() * self.x_t0d1() * self.x_t0d2()
            + 3.0 * self.site_x2::<C>(x) * self.x_t1d2() * self.x_t0d1()
            + 3.0 * self.site_x2::<C>(x) * self.x_t1d1() * self.x_t0d2()
            + self.site_x2::<C>(x) * self.x_t1d0() * self.x_t0d3()
            + self.site_x1::<C>(x) * self.x_t1d3()
    }
    fn site_t2d0<const C: i32>(&mut self, x: f64) -> f64 {
        self.site_x2::<C>(x) * self.x_t1d0().powi(2) + self.site_x1::<C>(x) * self.x_t2d0()
    }
}
impl PcSaftGlyPure {
    fn site_x1<const C: i32>(&self, x: f64) -> f64 {
        (1.0 / x - 0.5) * C as f64
    }
    fn site_x2<const C: i32>(&self, x: f64) -> f64 {
        -1.0 / x.powi(2) * C.pow(2) as f64
    }
    fn site_x3<const C: i32>(&self, x: f64) -> f64 {
        2.0 / x.powi(3) * C.pow(3) as f64
    }
    fn site_x4<const C: i32>(&self, x: f64) -> f64 {
        -6.0 / x.powi(4) * C.pow(4) as f64
    }
}
impl PcSaftGlyPure {
    fn xt1(&mut self) -> f64 {
        if self.xt1.0 != self.XA {
            self.xt1 = (
                self.XA,
                match self.assoc_type {
                    AssocType::Type2B => self.XA.powi(3) / (self.XA - 2.0),
                    AssocType::Type3B => {
                        (self.XA * (2.0 * self.XA - 1.0)).powi(2)
                            / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0)
                    }
                    AssocType::Type0 => 0.0,
                },
            )
        }
        self.xt1.1
    }
    fn xt2(&mut self) -> f64 {
        if self.xt2.0 != self.XA {
            self.xt2 = (
                self.XA,
                2.0 * match self.assoc_type {
                    AssocType::Type2B => {
                        self.XA.powi(5) / (self.XA - 2.0).powi(3) * (self.XA - 3.0)
                    }
                    AssocType::Type3B => {
                        (self.XA * (2.0 * self.XA - 1.0)).powi(3)
                            / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(3)
                            * (4.0 * self.XA.powi(3) - 12.0 * self.XA.powi(2) + 6.0 * self.XA - 1.0)
                    }
                    AssocType::Type0 => 0.0,
                },
            )
        }
        self.xt2.1
    }
    fn xt3(&mut self) -> f64 {
        if self.xt3.0 != self.XA {
            self.xt3 = (
                self.XA,
                6.0 * match self.assoc_type {
                    AssocType::Type2B => {
                        self.XA.powi(7) / (self.XA - 2.0).powi(5)
                            * (self.XA.powi(2) - 6.0 * self.XA + 10.0)
                    }
                    AssocType::Type3B => {
                        (self.XA * (2.0 * self.XA - 1.0)).powi(4)
                            / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(5)
                            * (16.0 * self.XA.powi(6) - 96.0 * self.XA.powi(5)
                                + (200.0 * self.XA.powi(4) - 160.0 * self.XA.powi(3))
                                + (62.0 * self.XA.powi(2) - 12.0 * self.XA + 1.0))
                    }
                    AssocType::Type0 => 0.0,
                },
            )
        }
        self.xt3.1
    }
    fn xt4(&mut self) -> f64 {
        if self.xt4.0 != self.XA {
            self.xt4 = (
                self.XA,
                24.0 * match self.assoc_type {
                    AssocType::Type2B => {
                        self.XA.powi(9) / (self.XA - 2.0).powi(7)
                            * (self.XA.powi(3) - 9.0 * self.XA.powi(2) + 29.0 * self.XA - 35.0)
                    }
                    AssocType::Type3B => {
                        (self.XA * (2.0 * self.XA - 1.0)).powi(5)
                            / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(7)
                            * (64.0 * self.XA.powi(9) - 576.0 * self.XA.powi(8)
                                + (2080.0 * self.XA.powi(7) - 3808.0 * self.XA.powi(6))
                                + (3696.0 * self.XA.powi(5) - 2084.0 * self.XA.powi(4))
                                + (716.0 * self.XA.powi(3) - 150.0 * self.XA.powi(2))
                                + (18.0 * self.XA - 1.0))
                    }
                    AssocType::Type0 => 0.0,
                },
            )
        }
        self.xt4.1
    }
}
impl PcSaftGlyPure {
    fn x_t0d1(&mut self) -> f64 {
        if self.x_t0d1.0 != self.XA {
            self.x_t0d1 = (self.XA, self.xt1() * self.t_t0d1())
        }
        self.x_t0d1.1
    }
    fn x_t0d2(&mut self) -> f64 {
        if self.x_t0d2.0 != self.XA {
            self.x_t0d2 = (
                self.XA,
                self.xt2() * self.t_t0d1().powi(2) + self.xt1() * self.t_t0d2(),
            )
        }
        self.x_t0d2.1
    }
    fn x_t0d3(&mut self) -> f64 {
        if self.x_t0d3.0 != self.XA {
            self.x_t0d3 = (
                self.XA,
                self.xt3() * self.t_t0d1().powi(3)
                    + 3.0 * self.xt2() * self.t_t0d1() * self.t_t0d2()
                    + self.xt1() * self.t_t0d3(),
            )
        }
        self.x_t0d3.1
    }
    fn x_t0d4(&mut self) -> f64 {
        if self.x_t0d4.0 != self.XA {
            self.x_t0d4 = (
                self.XA,
                self.xt4() * self.t_t0d1().powi(4)
                    + 6.0 * self.xt3() * self.t_t0d1().powi(2) * self.t_t0d2()
                    + 3.0 * self.xt2() * self.t_t0d2().powi(2)
                    + 4.0 * self.xt2() * self.t_t0d1() * self.t_t0d3()
                    + self.xt1() * self.t_t0d4(),
            )
        }
        self.x_t0d4.1
    }
    fn x_t1d0(&mut self) -> f64 {
        if self.x_t1d0.0 != self.XA {
            self.x_t1d0 = (self.XA, self.xt1() * self.t_t1d0())
        }
        self.x_t1d0.1
    }
    fn x_t1d1(&mut self) -> f64 {
        if self.x_t1d1.0 != self.XA {
            self.x_t1d1 = (
                self.XA,
                self.xt2() * self.t_t1d0() * self.t_t0d1() + self.xt1() * self.t_t1d1(),
            )
        }
        self.x_t1d1.1
    }
    fn x_t1d2(&mut self) -> f64 {
        if self.x_t1d2.0 != self.XA {
            self.x_t1d2 = (
                self.XA,
                self.xt3() * self.t_t1d0() * self.t_t0d1().powi(2)
                    + 2.0 * self.xt2() * self.t_t1d1() * self.t_t0d1()
                    + self.xt2() * self.t_t1d0() * self.t_t0d2()
                    + self.xt1() * self.t_t1d2(),
            )
        }
        self.x_t1d2.1
    }
    fn x_t1d3(&mut self) -> f64 {
        if self.x_t1d3.0 != self.XA {
            self.x_t1d3 = (
                self.XA,
                self.xt4() * self.t_t1d0() * self.t_t0d1().powi(3)
                    + 3.0 * self.xt3() * self.t_t1d1() * self.t_t0d1().powi(2)
                    + 3.0 * self.xt3() * self.t_t1d0() * self.t_t0d1() * self.t_t0d2()
                    + 3.0 * self.xt2() * self.t_t1d2() * self.t_t0d1()
                    + 3.0 * self.xt2() * self.t_t1d1() * self.t_t0d2()
                    + self.xt2() * self.t_t1d0() * self.t_t0d3()
                    + self.xt1() * self.t_t1d3(),
            )
        }
        self.x_t1d3.1
    }
    fn x_t2d0(&mut self) -> f64 {
        if self.x_t2d0.0 != self.XA {
            self.x_t2d0 = (
                self.XA,
                self.xt2() * self.t_t1d0().powi(2) + self.xt1() * self.t_t2d0(),
            )
        }
        self.x_t2d0.1
    }
}
impl PcSaftGlyPure {
    fn t_t0d0(&mut self) -> f64 {
        if self.eta0 != self.t_t0d0.0 {
            self.t_t0d0 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * self.gii_gly_t0d0(),
            )
        }
        self.t_t0d0.1
    }
    fn t_t0d1(&mut self) -> f64 {
        if self.eta0 != self.t_t0d1.0 {
            self.t_t0d1 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d1() + self.gii_gly_t0d0()),
            )
        }
        self.t_t0d1.1
    }
    fn t_t0d2(&mut self) -> f64 {
        if self.eta0 != self.t_t0d2.0 {
            self.t_t0d2 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d2() + 2.0 * self.gii_gly_t0d1()),
            )
        }
        self.t_t0d2.1
    }
    fn t_t0d3(&mut self) -> f64 {
        if self.eta0 != self.t_t0d3.0 {
            self.t_t0d3 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d3() + 3.0 * self.gii_gly_t0d2()),
            )
        }
        self.t_t0d3.1
    }
    fn t_t0d4(&mut self) -> f64 {
        if self.eta0 != self.t_t0d4.0 {
            self.t_t0d4 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d4() + 4.0 * self.gii_gly_t0d3()),
            )
        }
        self.t_t0d4.1
    }
    fn t_t1d0(&mut self) -> f64 {
        if self.eta0 != self.t_t1d0.0 {
            self.t_t1d0 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * ((self.epsilon_AB_temp.exp() - 1.0) * self.gii_gly_t1d0()
                        - self.epsilon_AB_temp.exp() * self.epsilon_AB_temp * self.gii_gly_t0d0()),
            )
        }
        self.t_t1d0.1
    }
    fn t_t1d1(&mut self) -> f64 {
        if self.eta0 != self.t_t1d1.0 {
            self.t_t1d1 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii_gly_t1d1() + self.gii_gly_t1d0())
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii_gly_t0d1() + self.gii_gly_t0d0())),
            )
        }
        self.t_t1d1.1
    }
    fn t_t1d2(&mut self) -> f64 {
        if self.eta0 != self.t_t1d2.0 {
            self.t_t1d2 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii_gly_t1d2() + 2.0 * self.gii_gly_t1d1())
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii_gly_t0d2() + 2.0 * self.gii_gly_t0d1())),
            )
        }
        self.t_t1d2.1
    }
    fn t_t1d3(&mut self) -> f64 {
        if self.eta0 != self.t_t1d3.0 {
            self.t_t1d3 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii_gly_t1d3() + 3.0 * self.gii_gly_t1d2())
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii_gly_t0d3() + 3.0 * self.gii_gly_t0d2())),
            )
        }
        self.t_t1d3.1
    }
    fn t_t2d0(&mut self) -> f64 {
        if self.eta0 != self.t_t2d0.0 {
            self.t_t2d0 = (
                self.eta0,
                (self.rho_num * self.kappa_AB_sigma3)
                    * (self.gii_gly_t2d0() * (self.epsilon_AB_temp.exp() - 1.0)
                        - 2.0
                            * self.gii_gly_t1d0()
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                        + self.gii_gly_t0d0()
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * (self.epsilon_AB_temp + 2.0)),
            )
        }
        self.t_t2d0.1
    }
}
impl PcSaftGlyPure {
    fn gii_gly_t0d0(&mut self) -> f64 {
        if self.eta0 != self.g_t0d0.0 {
            self.g_t0d0 = (
                self.eta0,
                (1.0 - self.eta0).powi(3).recip()
                    * (self.c0 + self.c1 * self.eta0 + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t0d0.1
    }
    fn gii_gly_t0d1(&mut self) -> f64 {
        if self.eta0 != self.g_t0d1.0 {
            self.g_t0d1 = (
                self.eta0,
                self.eta0 / (1.0 - self.eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t0d1.1
    }
    fn gii_gly_t0d2(&mut self) -> f64 {
        if self.eta0 != self.g_t0d2.0 {
            self.g_t0d2 = (
                self.eta0,
                2.0 * self.eta0.powi(2) / (1.0 - self.eta0).powi(5)
                    * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                        + (3.0 * self.c1 + 4.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t0d2.1
    }
    fn gii_gly_t0d3(&mut self) -> f64 {
        if self.eta0 != self.g_t0d3.0 {
            self.g_t0d3 = (
                self.eta0,
                6.0 * self.eta0.powi(3) / (1.0 - self.eta0).powi(6)
                    * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                        + (4.0 * self.c1 + 6.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t0d3.1
    }
    fn gii_gly_t0d4(&mut self) -> f64 {
        if self.eta0 != self.g_t0d4.0 {
            self.g_t0d4 = (
                self.eta0,
                24.0 * self.eta0.powi(4) / (1.0 - self.eta0).powi(7)
                    * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                        + (5.0 * self.c1 + 8.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t0d4.1
    }
    fn gii_gly_t1d0(&mut self) -> f64 {
        if self.eta0 != self.g_t1d0.0 {
            self.g_t1d0 = (
                self.eta0,
                self.eta1 / (1.0 - self.eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t1d0.1
    }
    fn gii_gly_t1d1(&mut self) -> f64 {
        if self.eta0 != self.g_t1d1.0 {
            self.g_t1d1 = (
                self.eta0,
                self.eta1 / (1.0 - self.eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2))
                    + 2.0 * self.eta1 * self.eta0 / (1.0 - self.eta0).powi(5)
                        * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                            + (3.0 * self.c1 + 4.0 * self.c2) * self.eta0
                            + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t1d1.1
    }
    fn gii_gly_t1d2(&mut self) -> f64 {
        if self.eta0 != self.g_t1d2.0 {
            self.g_t1d2 = (
                self.eta0,
                4.0 * self.eta1 * self.eta0 / (1.0 - self.eta0).powi(5)
                    * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                        + (3.0 * self.c1 + 4.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2))
                    + 6.0 * self.eta1 * self.eta0.powi(2) / (1.0 - self.eta0).powi(6)
                        * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                            + (4.0 * self.c1 + 6.0 * self.c2) * self.eta0
                            + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t1d2.1
    }
    fn gii_gly_t1d3(&mut self) -> f64 {
        if self.eta0 != self.g_t1d3.0 {
            self.g_t1d3 = (
                self.eta0,
                18.0 * self.eta1 * self.eta0.powi(2) / (1.0 - self.eta0).powi(6)
                    * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                        + (4.0 * self.c1 + 6.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2))
                    + 24.0 * self.eta1 * self.eta0.powi(3) / (1.0 - self.eta0).powi(7)
                        * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                            + (5.0 * self.c1 + 8.0 * self.c2) * self.eta0
                            + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t1d3.1
    }
    fn gii_gly_t2d0(&mut self) -> f64 {
        if self.eta0 != self.g_t2d0.0 {
            self.g_t2d0 = (
                self.eta0,
                self.eta2 / (1.0 - self.eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * self.eta0
                        + self.c2 * self.eta0.powi(2))
                    + 2.0 * self.eta1.powi(2) / (1.0 - self.eta0).powi(5)
                        * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                            + (3.0 * self.c1 + 4.0 * self.c2) * self.eta0
                            + self.c2 * self.eta0.powi(2)),
            )
        }
        self.g_t2d0.1
    }
}
