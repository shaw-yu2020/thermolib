use super::PcSaftErr;
use super::{AssocGlyTerm, DispTerm, GiiTerm, HsTerm};
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
/// methanol.set_2B_assoc_term(0.035176, 2899.5, 1.0, 1.0, 1.0); // kappa_AB epsilon_AB
/// methanol.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(methanol.rho().unwrap().round(), 24676.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)] // For pyclass hhh
pub struct PcSaftGlyPure {
    m: f64,
    sigma3: f64,
    epsilon: f64,
    hs: HsTerm,                  // HsTerm
    gii: GiiTerm,                // GiiTerm
    disp: DispTerm,              // DispTerm
    assoc: Option<AssocGlyTerm>, // AssocTerm
    // state
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
impl PcSaftGlyPure {
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
            // state
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
impl PcSaftGlyPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_1_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) {
        self.assoc = Some(AssocGlyTerm::new_1_term(
            kappa_AB * self.sigma3,
            epsilon_AB,
            c0,
            c1,
            c2,
        ));
    }
    pub fn set_2B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) {
        self.assoc = Some(AssocGlyTerm::new_2B_term(
            kappa_AB * self.sigma3,
            epsilon_AB,
            c0,
            c1,
            c2,
        ))
    }
    pub fn set_3B_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) {
        self.assoc = Some(AssocGlyTerm::new_3B_term(
            kappa_AB * self.sigma3,
            epsilon_AB,
            c0,
            c1,
            c2,
        ))
    }
    pub fn set_4C_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) {
        self.assoc = Some(AssocGlyTerm::new_4C_term(
            kappa_AB * self.sigma3,
            epsilon_AB,
            c0,
            c1,
            c2,
        ))
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
                (FRAC_PI_6 * self.m * self.sigma3)
                    * (3.0
                        * (1.0 - 0.12 * (-epsilon_temp_plus).exp()).powi(2)
                        * (-0.12 * (-epsilon_temp_plus).exp() * epsilon_temp_plus)),
            );
        }
        self.eta1_coef.1
    }
    fn eta2_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta2_coef.0 {
            let epsilon_temp_plus = 3.0 * self.epsilon / temp;
            self.eta2_coef = (
                temp,
                (FRAC_PI_6 * self.m * self.sigma3)
                    * (6.0
                        * (1.0 - 0.12 * (-epsilon_temp_plus).exp())
                        * (-0.12 * (-epsilon_temp_plus).exp() * epsilon_temp_plus).powi(2)
                        + 3.0
                            * (1.0 - 0.12 * (-epsilon_temp_plus).exp()).powi(2)
                            * (-0.12
                                * (-epsilon_temp_plus).exp()
                                * epsilon_temp_plus
                                * (epsilon_temp_plus - 2.0))),
            );
        }
        self.eta2_coef.1
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta0 = self.eta0_coef(temp) * rho_num;
        self.m * self.hs.t0d0(eta0) - (self.m - 1.0) * self.gii.lngii_t0d0(eta0)
            + self.disp.t0d0(temp, rho_num, eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d0(temp, rho_num, eta0))
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta0 = self.eta0_coef(temp) * rho_num;
        self.m * self.hs.t0d1(eta0) - (self.m - 1.0) * self.gii.lngii_t0d1(eta0)
            + self.disp.t0d1(temp, rho_num, eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d1(temp, rho_num, eta0))
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta0 = self.eta0_coef(temp) * rho_num;
        self.m * self.hs.t0d2(eta0) - (self.m - 1.0) * self.gii.lngii_t0d2(eta0)
            + self.disp.t0d2(temp, rho_num, eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d2(temp, rho_num, eta0))
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta0 = self.eta0_coef(temp) * rho_num;
        self.m * self.hs.t0d3(eta0) - (self.m - 1.0) * self.gii.lngii_t0d3(eta0)
            + self.disp.t0d3(temp, rho_num, eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d3(temp, rho_num, eta0))
    }
    fn r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta0 = self.eta0_coef(temp) * rho_num;
        self.m * self.hs.t0d4(eta0) - (self.m - 1.0) * self.gii.lngii_t0d4(eta0)
            + self.disp.t0d4(temp, rho_num, eta0)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t0d4(temp, rho_num, eta0))
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta0, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.m * self.hs.t1d0(eta0, eta1) - (self.m - 1.0) * self.gii.lngii_t1d0(eta0, eta1)
            + self.disp.t1d0(temp, rho_num, eta0, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d0(temp, rho_num, eta0, eta1))
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta0, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.m * self.hs.t1d1(eta0, eta1) - (self.m - 1.0) * self.gii.lngii_t1d1(eta0, eta1)
            + self.disp.t1d1(temp, rho_num, eta0, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d1(temp, rho_num, eta0, eta1))
    }
    fn r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta0, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.m * self.hs.t1d2(eta0, eta1) - (self.m - 1.0) * self.gii.lngii_t1d2(eta0, eta1)
            + self.disp.t1d2(temp, rho_num, eta0, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d2(temp, rho_num, eta0, eta1))
    }
    fn r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta0, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.m * self.hs.t1d3(eta0, eta1) - (self.m - 1.0) * self.gii.lngii_t1d3(eta0, eta1)
            + self.disp.t1d3(temp, rho_num, eta0, eta1)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t1d3(temp, rho_num, eta0, eta1))
    }
    fn r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta0, eta1, eta2) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
            self.eta2_coef(temp) * rho_num,
        );
        self.m * self.hs.t2d0(eta0, eta1, eta2)
            - (self.m - 1.0) * self.gii.lngii_t2d0(eta0, eta1, eta2)
            + self.disp.t2d0(temp, rho_num, eta0, eta1, eta2)
            + self
                .assoc
                .as_mut()
                .map_or(0.0, |a| a.t2d0(temp, rho_num, eta0, eta1, eta2))
    }
}
