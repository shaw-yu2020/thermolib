use super::PcSaftErr;
use super::{AssocTerm, DispTerm, GiiTerm, HsTerm, PddTerm, PqqTerm};
use super::{FRAC_NA_1E30, FRAC_RE30_NA, R};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::FRAC_PI_6;
/// PC-SAFT EOS
/// ```
/// use thermolib::PcSaftPure;
/// let (m, sigma, epsilon) = (2.8611, 2.6826, 205.35); // SO2
/// let mut fluid = PcSaftPure::new_fluid(m, sigma, epsilon);
/// fluid.c_flash().unwrap();
/// assert_eq!(fluid.T().unwrap().round(), 438.0);
/// assert_eq!(fluid.p().unwrap().round(), 9099261.0);
/// assert_eq!(fluid.rho().unwrap().round(), 8079.0);
/// fluid.t_flash(298.15).unwrap();
/// assert_eq!(fluid.p_s().unwrap().round(), 396865.0);
/// assert_eq!(fluid.rho_v().unwrap().round(), 169.0);
/// assert_eq!(fluid.rho_l().unwrap().round(), 21115.0);
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 41.0);
/// let mut fluid = PcSaftPure::new_fluid(1.0656, 3.0007, 366.51); // H2O
/// fluid.set_2B_assoc_term(0.034868, 2500.7); // kappa_AB epsilon_AB
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 51179.0);
/// let mut fluid = PcSaftPure::new_fluid(1.5131, 3.1869, 163.33); // CO2
/// fluid.set_QQ_polar_term(4.4); // |Q|(DA)
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 41.0);
/// let mut fluid = PcSaftPure::new_fluid(2.7447, 3.2742, 232.99); // Acetone
/// fluid.set_DD_polar_term(2.88); // |Q|(DA)
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 13276.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)] // For pyclass hhh
pub struct PcSaftPure {
    m: f64,
    sigma3: f64,
    epsilon: f64,
    hs: HsTerm,               // HsTerm
    gii: GiiTerm,             // GiiTerm
    disp: DispTerm,           // DispTerm
    assoc: Option<AssocTerm>, // AssocTerm
    qq: Option<PqqTerm>,      // PqqTerm
    dd: Option<PddTerm>,      // PddTerm
    // state
    eta0: f64,
    eta1: f64,
    eta2: f64,
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
}
impl PcSaftPure {
    pub fn new_fluid(m: f64, sigma: f64, epsilon: f64) -> Self {
        let sigma3 = sigma.powi(3);
        Self {
            m,
            sigma3,
            epsilon,
            hs: HsTerm::new(),                       // HsTerm
            gii: GiiTerm::new(),                     // GiiTerm
            disp: DispTerm::new(m, sigma3, epsilon), // DispTerm
            assoc: None,                             // AssocTerm
            dd: None,                                // PqqTerm
            qq: None,                                // PddTerm
            // state
            eta0: 0.0,
            eta1: 0.0,
            eta2: 0.0,
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
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)] // For pymethods hhh
impl PcSaftPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_1_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_1_term(kappa_AB * self.sigma3, epsilon_AB));
    }
    pub fn set_2B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_2B_term(kappa_AB * self.sigma3, epsilon_AB))
    }
    pub fn set_3B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_3B_term(kappa_AB * self.sigma3, epsilon_AB))
    }
    pub fn set_4C_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocTerm::new_4C_term(kappa_AB * self.sigma3, epsilon_AB))
    }
    pub fn set_QQ_polar_term(&mut self, Q: f64) {
        self.qq = Some(PqqTerm::new(self.m, self.sigma3, self.epsilon, Q))
    }
    pub fn set_DD_polar_term(&mut self, u: f64) {
        self.dd = Some(PddTerm::new(self.m, self.sigma3, self.epsilon, u))
    }
}
fn_vec!(PcSaftPure);
fn_flash!(PcSaftPure);
fn_c_flash!(PcSaftPure);
fn_single_prop!(PcSaftPure);
fn_double_prop!(PcSaftPure);
fn_virial_prop!(PcSaftPure);
fn_aly_lee_cp0!(PcSaftPure);
impl PcSaftPure {
    fn_calc_prop!();
    fn_check_derivatives!();
}
impl PcSaftPure {
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
        }
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d0(self.eta0) - (self.m - 1.0) * self.gii.lngii_t0d0(self.eta0)
            + self.disp.t0d0(temp, rho_num, self.eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d0(temp, rho_num, self.eta0))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t0d0(temp, rho_num, self.eta0))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t0d0(temp, rho_num, self.eta0))
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d1(self.eta0) - (self.m - 1.0) * self.gii.lngii_t0d1(self.eta0)
            + self.disp.t0d1(temp, rho_num, self.eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d1(temp, rho_num, self.eta0))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t0d1(temp, rho_num, self.eta0))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t0d1(temp, rho_num, self.eta0))
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d2(self.eta0) - (self.m - 1.0) * self.gii.lngii_t0d2(self.eta0)
            + self.disp.t0d2(temp, rho_num, self.eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d2(temp, rho_num, self.eta0))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t0d2(temp, rho_num, self.eta0))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t0d2(temp, rho_num, self.eta0))
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d3(self.eta0) - (self.m - 1.0) * self.gii.lngii_t0d3(self.eta0)
            + self.disp.t0d3(temp, rho_num, self.eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d3(temp, rho_num, self.eta0))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t0d3(temp, rho_num, self.eta0))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t0d3(temp, rho_num, self.eta0))
    }
    fn r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t0d4(self.eta0) - (self.m - 1.0) * self.gii.lngii_t0d4(self.eta0)
            + self.disp.t0d4(temp, rho_num, self.eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d4(temp, rho_num, self.eta0))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t0d4(temp, rho_num, self.eta0))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t0d4(temp, rho_num, self.eta0))
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d0(self.eta0, self.eta1)
            - (self.m - 1.0) * self.gii.lngii_t1d0(self.eta0, self.eta1)
            + self.disp.t1d0(temp, rho_num, self.eta0, self.eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d0(temp, rho_num, self.eta0, self.eta1))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t1d0(temp, rho_num, self.eta0, self.eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t1d0(temp, rho_num, self.eta0, self.eta1))
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d1(self.eta0, self.eta1)
            - (self.m - 1.0) * self.gii.lngii_t1d1(self.eta0, self.eta1)
            + self.disp.t1d1(temp, rho_num, self.eta0, self.eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d1(temp, rho_num, self.eta0, self.eta1))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t1d1(temp, rho_num, self.eta0, self.eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t1d1(temp, rho_num, self.eta0, self.eta1))
    }
    fn r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d2(self.eta0, self.eta1)
            - (self.m - 1.0) * self.gii.lngii_t1d2(self.eta0, self.eta1)
            + self.disp.t1d2(temp, rho_num, self.eta0, self.eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d2(temp, rho_num, self.eta0, self.eta1))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t1d2(temp, rho_num, self.eta0, self.eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t1d2(temp, rho_num, self.eta0, self.eta1))
    }
    fn r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t1d3(self.eta0, self.eta1)
            - (self.m - 1.0) * self.gii.lngii_t1d3(self.eta0, self.eta1)
            + self.disp.t1d3(temp, rho_num, self.eta0, self.eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d3(temp, rho_num, self.eta0, self.eta1))
            + self
                .qq
                .as_mut()
                .map_or(0.0, |a| a.t1d3(temp, rho_num, self.eta0, self.eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |a| a.t1d3(temp, rho_num, self.eta0, self.eta1))
    }
    fn r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.eta_flash(temp, rho_num);
        self.m * self.hs.t2d0(self.eta0, self.eta1, self.eta2)
            - (self.m - 1.0) * self.gii.lngii_t2d0(self.eta0, self.eta1, self.eta2)
            + self
                .disp
                .t2d0(temp, rho_num, self.eta0, self.eta1, self.eta2)
            + self.assoc.as_mut().map_or(0.0, |a| {
                a.t2d0(temp, rho_num, self.eta0, self.eta1, self.eta2)
            })
            + self.qq.as_mut().map_or(0.0, |a| {
                a.t2d0(temp, rho_num, self.eta0, self.eta1, self.eta2)
            })
            + self.dd.as_mut().map_or(0.0, |a| {
                a.t2d0(temp, rho_num, self.eta0, self.eta1, self.eta2)
            })
    }
}
