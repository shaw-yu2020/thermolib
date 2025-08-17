use super::PcSaftErr;
use super::{AssocPure, DispTerm, GiiPure, HsPure, PolarTerm};
use super::{FRAC_NA_1E30, FRAC_RE30_NA, R};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_6};
/// PC-SAFT EOS :: PcSaftPure
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
/// let mut fluid = PcSaftPure::new_fluid(1.0656, 3.0007, 366.51); // H2O_2B
/// fluid.set_2B_assoc_term(0.034868, 2500.7); // kappa_AB epsilon_AB
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 51179.0);
/// let mut fluid = PcSaftPure::new_fluid(1.0656, 3.0007, 366.51); // H2O_3B
/// fluid.set_3B_assoc_term(0.034868, 2500.7); // kappa_AB epsilon_AB
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 51371.0);
/// let mut fluid = PcSaftPure::new_fluid(1.0656, 3.0007, 366.51); // H2O_4C
/// fluid.set_4C_assoc_term(0.034868, 2500.7); // kappa_AB epsilon_AB
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 55670.0);
/// let mut fluid = PcSaftPure::new_fluid(1.5131, 3.1869, 163.33); // CO2_QQ
/// fluid.set_QQ_polar_term(4.4); // |Q|(DA)
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 41.0);
/// let mut fluid = PcSaftPure::new_fluid(2.7447, 3.2742, 232.99); // ACETONE_DD
/// fluid.set_DD_polar_term(2.88); // |u|(D)
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 13337.0);
/// let mut fluid = PcSaftPure::new_fluid(1.5255, 3.23, 188.9); // METHANOL_2B
/// fluid.set_2B_assoc_term(0.035176, 2899.5); // kappa_AB epsilon_AB
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 24676.0);
/// let mut fluid = PcSaftPure::new_fluid(1.5255, 3.23, 188.9); // METHANOL_3B
/// fluid.set_3B_assoc_term(0.035176, 2899.5); // kappa_AB epsilon_AB
/// fluid.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 24836.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)] // For pyclass hhh
pub struct PcSaftPure {
    m: f64,
    sigma3: f64,
    epsilon: f64,
    hs: HsPure,               // HsPure
    gii: GiiPure,             // GiiPure
    disp: DispTerm,           // DispTerm
    assoc: Option<AssocPure>, // AssocPure
    polar: Option<PolarTerm>, // PolarTerm
    temp: f64,
    rho_num: f64,
    rhov_num: f64,
    rhol_num: f64,
    eta0_coef: (f64, f64),
    eta1_coef: (f64, f64),
    eta2_coef: (f64, f64),
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
            hs: HsPure::new(m),         // HsPure
            gii: GiiPure::new(m - 1.0), // GiiPure
            disp: DispTerm::new(
                m,
                m.powi(2) * epsilon * sigma3,
                m.powi(2) * epsilon.powi(2) * sigma3,
            ), // DispTerm
            assoc: None,                // AssocPure
            polar: None,                // PolarTerm
            temp: 1.0,
            rho_num: 1.0,
            rhov_num: 0.0,
            rhol_num: 0.0,
            eta0_coef: (0.0, 0.0),
            eta1_coef: (0.0, 0.0),
            eta2_coef: (0.0, 0.0),
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
        self.assoc = Some(AssocPure::new_1_term(kappa_AB * self.sigma3, epsilon_AB));
    }
    pub fn set_2B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocPure::new_2B_term(kappa_AB * self.sigma3, epsilon_AB))
    }
    pub fn set_3B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocPure::new_3B_term(kappa_AB * self.sigma3, epsilon_AB))
    }
    pub fn set_4C_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc = Some(AssocPure::new_4C_term(kappa_AB * self.sigma3, epsilon_AB))
    }
    pub fn set_QQ_polar_term(&mut self, Q: f64) {
        self.polar = Some(PolarTerm::new(
            &[1.0],
            &[self.m],
            &[self.sigma3.cbrt()],
            &[self.epsilon],
            &[Q],
            &[4],
        ));
    }
    pub fn set_DD_polar_term(&mut self, u: f64) {
        self.polar = Some(PolarTerm::new(
            &[1.0],
            &[self.m],
            &[self.sigma3.cbrt()],
            &[self.epsilon],
            &[u],
            &[2],
        ))
    }
}
fn_vec!(PcSaftPure);
fn_c_flash!(PcSaftPure);
fn_t_flash!(PcSaftPure);
fn_tp_flash!(PcSaftPure);
fn_single_prop!(PcSaftPure);
fn_double_prop!(PcSaftPure);
fn_virial_prop!(PcSaftPure);
fn_aly_lee_cp0!(PcSaftPure);
impl PcSaftPure {
    fn_calc_prop!();
    fn_check_derivatives!();
}
impl PcSaftPure {
    fn eta0_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta0_coef.0 {
            self.eta0_coef = (
                temp,
                (FRAC_PI_6 * self.m * self.sigma3)
                    * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp()).powi(3),
            );
        }
        self.eta0_coef.1
    }
    fn eta1_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta1_coef.0 {
            let epsilon_temp_plus = 3.0 * self.epsilon / temp;
            self.eta1_coef = (
                temp,
                (FRAC_PI_2 * self.m * self.sigma3)
                    * (1.0 - 0.12 * (-epsilon_temp_plus).exp()).powi(2)
                    * (-0.12 * (-epsilon_temp_plus).exp() * epsilon_temp_plus),
            );
        }
        self.eta1_coef.1
    }
    fn eta2_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta2_coef.0 {
            let epsilon_temp_plus = 3.0 * self.epsilon / temp;
            self.eta2_coef = (
                temp,
                (FRAC_PI_2 * self.m * self.sigma3)
                    * (2.0
                        * (1.0 - 0.12 * (-epsilon_temp_plus).exp())
                        * (-0.12 * (-epsilon_temp_plus).exp() * epsilon_temp_plus).powi(2)
                        + (1.0 - 0.12 * (-epsilon_temp_plus).exp()).powi(2)
                            * (-0.12
                                * (-epsilon_temp_plus).exp()
                                * epsilon_temp_plus
                                * (epsilon_temp_plus - 2.0))),
            );
        }
        self.eta2_coef.1
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d0(eta)
            + self.gii.lngii_t0d0(eta)
            + self.disp.t0d0(temp, rho_num, eta)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d0(temp, rho_num, eta))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d0(temp, rho_num, eta))
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d1(eta)
            + self.gii.lngii_t0d1(eta)
            + self.disp.t0d1(temp, rho_num, eta)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d1(temp, rho_num, eta))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d1(temp, rho_num, eta))
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d2(eta)
            + self.gii.lngii_t0d2(eta)
            + self.disp.t0d2(temp, rho_num, eta)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d2(temp, rho_num, eta))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d2(temp, rho_num, eta))
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d3(eta)
            + self.gii.lngii_t0d3(eta)
            + self.disp.t0d3(temp, rho_num, eta)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d3(temp, rho_num, eta))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d3(temp, rho_num, eta))
    }
    fn r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d4(eta)
            + self.gii.lngii_t0d4(eta)
            + self.disp.t0d4(temp, rho_num, eta)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d4(temp, rho_num, eta))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t0d4(temp, rho_num, eta))
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d0(eta, eta1)
            + self.gii.lngii_t1d0(eta, eta1)
            + self.disp.t1d0(temp, rho_num, eta, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d0(temp, rho_num, eta, eta1))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t1d0(temp, rho_num, eta, eta1))
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d1(eta, eta1)
            + self.gii.lngii_t1d1(eta, eta1)
            + self.disp.t1d1(temp, rho_num, eta, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d1(temp, rho_num, eta, eta1))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t1d1(temp, rho_num, eta, eta1))
    }
    fn r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d2(eta, eta1)
            + self.gii.lngii_t1d2(eta, eta1)
            + self.disp.t1d2(temp, rho_num, eta, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d2(temp, rho_num, eta, eta1))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t1d2(temp, rho_num, eta, eta1))
    }
    fn r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d3(eta, eta1)
            + self.gii.lngii_t1d3(eta, eta1)
            + self.disp.t1d3(temp, rho_num, eta, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d3(temp, rho_num, eta, eta1))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t1d3(temp, rho_num, eta, eta1))
    }
    fn r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1, eta2) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
            self.eta2_coef(temp) * rho_num,
        );
        self.hs.t2d0(eta, eta1, eta2)
            + self.gii.lngii_t2d0(eta, eta1, eta2)
            + self.disp.t2d0(temp, rho_num, eta, eta1, eta2)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t2d0(temp, rho_num, eta, eta1, eta2))
            + self
                .polar
                .as_mut()
                .map_or(0.0, |p| p.t2d0(temp, rho_num, eta, eta1, eta2))
    }
}
