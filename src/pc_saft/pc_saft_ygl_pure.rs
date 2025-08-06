use super::{DispTerm, GiiPure, HsPure, PcSaftErr};
use super::{FRAC_NA_1E30, FRAC_RE30_NA, R};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_6};
/// PC-SAFT EOS
/// ```
/// use thermolib::PcSaftYglPure;
/// let (m, sigma, epsilon) = (2.8611, 2.6826, 205.35); // SO2
/// let mut fluid = PcSaftYglPure::new_fluid(m, sigma, epsilon);
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
/// let mut methanol = PcSaftYglPure::new_fluid(1.5255, 3.23, 188.9);
/// methanol.set_ygl_assoc_term(0.035176, 2899.5, 0.0); // kappa_AB epsilon_AB
/// methanol.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(methanol.rho().unwrap().round(), 24676.0);
/// methanol.set_ygl_assoc_term(0.035176, 2899.5, 1.0); // kappa_AB epsilon_AB
/// methanol.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(methanol.rho().unwrap().round(), 24836.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)] // For pyclass hhh
pub struct PcSaftYglPure {
    m: f64,
    sigma3: f64,
    epsilon: f64,
    hs: HsPure,     // HsPure
    gii: GiiPure,   // GiiPure
    disp: DispTerm, // DispTerm
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
    // association term
    no_assoc_term: bool,
    kappa_AB_sigma3_dens: f64,
    epsilon_AB_temp: f64,
    kappa_AB_sigma3: f64,
    epsilon_AB: f64,
    b_assoc: f64,
    XA: (f64, f64),     // (eta, XA)
    xt1: (f64, f64),    // cached variables
    xt2: (f64, f64),    // cached variables
    xt3: (f64, f64),    // cached variables
    xt4: (f64, f64),    // cached variables
    x_t0d1: (f64, f64), // cached variables
    x_t0d2: (f64, f64), // cached variables
    x_t0d3: (f64, f64), // cached variables
    x_t0d4: (f64, f64), // cached variables
    x_t1d0: (f64, f64), // cached variables
    x_t1d1: (f64, f64), // cached variables
    x_t1d2: (f64, f64), // cached variables
    x_t1d3: (f64, f64), // cached variables
    x_t2d0: (f64, f64), // cached variables
    t_t0d0: (f64, f64), // cached variables
    t_t0d1: (f64, f64), // cached variables
    t_t0d2: (f64, f64), // cached variables
    t_t0d3: (f64, f64), // cached variables
    t_t0d4: (f64, f64), // cached variables
    t_t1d0: (f64, f64), // cached variables
    t_t1d1: (f64, f64), // cached variables
    t_t1d2: (f64, f64), // cached variables
    t_t1d3: (f64, f64), // cached variables
    t_t2d0: (f64, f64), // cached variables
}
impl PcSaftYglPure {
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
            // association term
            no_assoc_term: true,
            kappa_AB_sigma3_dens: 0.0,
            epsilon_AB_temp: 0.0,
            kappa_AB_sigma3: 0.0,
            epsilon_AB: 0.0,
            b_assoc: 0.0,
            XA: (0.0, 1.0),
            xt1: (0.0, 0.0),    // cached variables
            xt2: (0.0, 0.0),    // cached variables
            xt3: (0.0, 0.0),    // cached variables
            xt4: (0.0, 0.0),    // cached variables
            x_t0d1: (0.0, 0.0), // cached variables
            x_t0d2: (0.0, 0.0), // cached variables
            x_t0d3: (0.0, 0.0), // cached variables
            x_t0d4: (0.0, 0.0), // cached variables
            x_t1d0: (0.0, 0.0), // cached variables
            x_t1d1: (0.0, 0.0), // cached variables
            x_t1d2: (0.0, 0.0), // cached variables
            x_t1d3: (0.0, 0.0), // cached variables
            x_t2d0: (0.0, 0.0), // cached variables
            t_t0d0: (0.0, 0.0), // cached variables
            t_t0d1: (0.0, 0.0), // cached variables
            t_t0d2: (0.0, 0.0), // cached variables
            t_t0d3: (0.0, 0.0), // cached variables
            t_t0d4: (0.0, 0.0), // cached variables
            t_t1d0: (0.0, 0.0), // cached variables
            t_t1d1: (0.0, 0.0), // cached variables
            t_t1d2: (0.0, 0.0), // cached variables
            t_t1d3: (0.0, 0.0), // cached variables
            t_t2d0: (0.0, 0.0), // cached variables
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)] // For pymethods hhh
impl PcSaftYglPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn py_new(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_ygl_assoc_term(&mut self, kappa_AB: f64, epsilon_AB: f64, b_assoc: f64) {
        self.no_assoc_term = false;
        self.kappa_AB_sigma3 = kappa_AB * self.sigma3;
        self.epsilon_AB = epsilon_AB;
        self.b_assoc = b_assoc;
    }
}
fn_vec!(PcSaftYglPure);
fn_c_flash!(PcSaftYglPure);
fn_t_flash!(PcSaftYglPure);
fn_tp_flash!(PcSaftYglPure);
fn_single_prop!(PcSaftYglPure);
fn_double_prop!(PcSaftYglPure);
fn_virial_prop!(PcSaftYglPure);
fn_aly_lee_cp0!(PcSaftYglPure);
impl PcSaftYglPure {
    fn_calc_prop!();
    fn_check_derivatives!();
}
impl PcSaftYglPure {
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
    #[allow(non_snake_case)]
    fn XA_flash(&mut self, temp: f64, rho_num: f64, eta: f64) {
        if eta != self.XA.0 {
            self.kappa_AB_sigma3_dens = self.kappa_AB_sigma3 * rho_num;
            self.epsilon_AB_temp = self.epsilon_AB / temp;
            self.XA = (
                eta,
                (-(1.0 - self.b_assoc * self.t_t0d0(eta))
                    + ((1.0 - self.b_assoc * self.t_t0d0(eta)).powi(2)
                        + 4.0 * (1.0 + self.b_assoc) * self.t_t0d0(eta))
                    .sqrt())
                    / (2.0 * (1.0 + self.b_assoc) * self.t_t0d0(eta)),
            );
        }
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d0(eta)
            + self.gii.lngii_t0d0(eta)
            + self.disp.t0d0(temp, rho_num, eta)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.XA.1.ln()
                    + (self.b_assoc * self.XA.1 + 1.0 - self.b_assoc).ln()
                    + ((1.0 + self.b_assoc) * self.XA.1 - self.b_assoc).ln()
                    - (1.0 + self.b_assoc) * self.XA.1
                    + (1.0 + self.b_assoc)
            }
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d1(eta)
            + self.gii.lngii_t0d1(eta)
            + self.disp.t0d1(temp, rho_num, eta)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t0d1(1.0, self.XA.1, eta)
                    + self.site_t0d1(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                    )
                    + self.site_t0d1(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                    )
            }
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d2(eta)
            + self.gii.lngii_t0d2(eta)
            + self.disp.t0d2(temp, rho_num, eta)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t0d2(1.0, self.XA.1, eta)
                    + self.site_t0d2(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                    )
                    + self.site_t0d2(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                    )
            }
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d3(eta)
            + self.gii.lngii_t0d3(eta)
            + self.disp.t0d3(temp, rho_num, eta)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t0d3(1.0, self.XA.1, eta)
                    + self.site_t0d3(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                    )
                    + self.site_t0d3(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                    )
            }
    }
    fn r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d4(eta)
            + self.gii.lngii_t0d4(eta)
            + self.disp.t0d4(temp, rho_num, eta)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t0d4(1.0, self.XA.1, eta)
                    + self.site_t0d4(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                    )
                    + self.site_t0d4(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                    )
            }
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d0(eta, eta1)
            + self.gii.lngii_t1d0(eta, eta1)
            + self.disp.t1d0(temp, rho_num, eta, eta1)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t1d0(1.0, self.XA.1, eta, eta1)
                    + self.site_t1d0(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                        eta1,
                    )
                    + self.site_t1d0(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                        eta1,
                    )
            }
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d1(eta, eta1)
            + self.gii.lngii_t1d1(eta, eta1)
            + self.disp.t1d1(temp, rho_num, eta, eta1)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t1d1(1.0, self.XA.1, eta, eta1)
                    + self.site_t1d1(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                        eta1,
                    )
                    + self.site_t1d1(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                        eta1,
                    )
            }
    }
    fn r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d2(eta, eta1)
            + self.gii.lngii_t1d2(eta, eta1)
            + self.disp.t1d2(temp, rho_num, eta, eta1)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t1d2(1.0, self.XA.1, eta, eta1)
                    + self.site_t1d2(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                        eta1,
                    )
                    + self.site_t1d2(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                        eta1,
                    )
            }
    }
    fn r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d3(eta, eta1)
            + self.gii.lngii_t1d3(eta, eta1)
            + self.disp.t1d3(temp, rho_num, eta, eta1)
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t1d3(1.0, self.XA.1, eta, eta1)
                    + self.site_t1d3(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                        eta1,
                    )
                    + self.site_t1d3(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                        eta1,
                    )
            }
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
            + if self.no_assoc_term {
                0.0
            } else {
                self.XA_flash(temp, rho_num, eta);
                self.site_t2d0(1.0, self.XA.1, eta, eta1, eta2)
                    + self.site_t2d0(
                        self.b_assoc,
                        self.b_assoc * self.XA.1 + 1.0 - self.b_assoc,
                        eta,
                        eta1,
                        eta2,
                    )
                    + self.site_t2d0(
                        1.0 + self.b_assoc,
                        (1.0 + self.b_assoc) * self.XA.1 - self.b_assoc,
                        eta,
                        eta1,
                        eta2,
                    )
            }
    }
}
impl PcSaftYglPure {
    fn site_t0d1(&mut self, c: f64, x: f64, eta: f64) -> f64 {
        c * self.site_x1(x) * self.x_t0d1(eta)
    }
    fn site_t0d2(&mut self, c: f64, x: f64, eta: f64) -> f64 {
        c.powi(2) * self.site_x2(x) * self.x_t0d1(eta).powi(2)
            + c * self.site_x1(x) * self.x_t0d2(eta)
    }
    fn site_t0d3(&mut self, c: f64, x: f64, eta: f64) -> f64 {
        c.powi(3) * self.site_x3(x) * self.x_t0d1(eta).powi(3)
            + 3.0 * c.powi(2) * self.site_x2(x) * self.x_t0d1(eta) * self.x_t0d2(eta)
            + c * self.site_x1(x) * self.x_t0d3(eta)
    }
    fn site_t0d4(&mut self, c: f64, x: f64, eta: f64) -> f64 {
        c.powi(4) * self.site_x4(x) * self.x_t0d1(eta).powi(4)
            + 6.0 * c.powi(3) * self.site_x3(x) * self.x_t0d1(eta).powi(2) * self.x_t0d2(eta)
            + 3.0 * c.powi(2) * self.site_x2(x) * self.x_t0d2(eta).powi(2)
            + 4.0 * c.powi(2) * self.site_x2(x) * self.x_t0d1(eta) * self.x_t0d3(eta)
            + c * self.site_x1(x) * self.x_t0d4(eta)
    }
    fn site_t1d0(&mut self, c: f64, x: f64, eta: f64, eta1: f64) -> f64 {
        c * self.site_x1(x) * self.x_t1d0(eta, eta1)
    }
    fn site_t1d1(&mut self, c: f64, x: f64, eta: f64, eta1: f64) -> f64 {
        c.powi(2) * self.site_x2(x) * self.x_t1d0(eta, eta1) * self.x_t0d1(eta)
            + c * self.site_x1(x) * self.x_t1d1(eta, eta1)
    }
    fn site_t1d2(&mut self, c: f64, x: f64, eta: f64, eta1: f64) -> f64 {
        c.powi(3) * self.site_x3(x) * self.x_t1d0(eta, eta1) * self.x_t0d1(eta).powi(2)
            + 2.0 * c.powi(2) * self.site_x2(x) * self.x_t1d1(eta, eta1) * self.x_t0d1(eta)
            + c.powi(2) * self.site_x2(x) * self.x_t1d0(eta, eta1) * self.x_t0d2(eta)
            + c * self.site_x1(x) * self.x_t1d2(eta, eta1)
    }
    fn site_t1d3(&mut self, c: f64, x: f64, eta: f64, eta1: f64) -> f64 {
        c.powi(4) * self.site_x4(x) * self.x_t1d0(eta, eta1) * self.x_t0d1(eta).powi(3)
            + 3.0 * c.powi(3) * self.site_x3(x) * self.x_t1d1(eta, eta1) * self.x_t0d1(eta).powi(2)
            + 3.0
                * c.powi(3)
                * self.site_x3(x)
                * self.x_t1d0(eta, eta1)
                * self.x_t0d1(eta)
                * self.x_t0d2(eta)
            + 3.0 * c.powi(2) * self.site_x2(x) * self.x_t1d2(eta, eta1) * self.x_t0d1(eta)
            + 3.0 * c.powi(2) * self.site_x2(x) * self.x_t1d1(eta, eta1) * self.x_t0d2(eta)
            + c.powi(2) * self.site_x2(x) * self.x_t1d0(eta, eta1) * self.x_t0d3(eta)
            + c * self.site_x1(x) * self.x_t1d3(eta, eta1)
    }
    fn site_t2d0(&mut self, c: f64, x: f64, eta: f64, eta1: f64, eta2: f64) -> f64 {
        c.powi(2) * self.site_x2(x) * self.x_t1d0(eta, eta1).powi(2)
            + c * self.site_x1(x) * self.x_t2d0(eta, eta1, eta2)
    }
}
impl PcSaftYglPure {
    fn site_x1(&self, x: f64) -> f64 {
        1.0 / x - 0.5
    }
    fn site_x2(&self, x: f64) -> f64 {
        -1.0 / x.powi(2)
    }
    fn site_x3(&self, x: f64) -> f64 {
        2.0 / x.powi(3)
    }
    fn site_x4(&self, x: f64) -> f64 {
        -6.0 / x.powi(4)
    }
}
impl PcSaftYglPure {
    fn x_t0d1(&mut self, eta: f64) -> f64 {
        if self.x_t0d1.0 != self.XA.1 {
            self.x_t0d1 = (self.XA.1, self.xt1() * self.t_t0d1(eta))
        }
        self.x_t0d1.1
    }
    fn x_t0d2(&mut self, eta: f64) -> f64 {
        if self.x_t0d2.0 != self.XA.1 {
            self.x_t0d2 = (
                self.XA.1,
                self.xt2() * self.t_t0d1(eta).powi(2) + self.xt1() * self.t_t0d2(eta),
            )
        }
        self.x_t0d2.1
    }
    fn x_t0d3(&mut self, eta: f64) -> f64 {
        if self.x_t0d3.0 != self.XA.1 {
            self.x_t0d3 = (
                self.XA.1,
                self.xt3() * self.t_t0d1(eta).powi(3)
                    + 3.0 * self.xt2() * self.t_t0d1(eta) * self.t_t0d2(eta)
                    + self.xt1() * self.t_t0d3(eta),
            )
        }
        self.x_t0d3.1
    }
    fn x_t0d4(&mut self, eta: f64) -> f64 {
        if self.x_t0d4.0 != self.XA.1 {
            self.x_t0d4 = (
                self.XA.1,
                self.xt4() * self.t_t0d1(eta).powi(4)
                    + 6.0 * self.xt3() * self.t_t0d1(eta).powi(2) * self.t_t0d2(eta)
                    + 3.0 * self.xt2() * self.t_t0d2(eta).powi(2)
                    + 4.0 * self.xt2() * self.t_t0d1(eta) * self.t_t0d3(eta)
                    + self.xt1() * self.t_t0d4(eta),
            )
        }
        self.x_t0d4.1
    }
    fn x_t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        if self.x_t1d0.0 != self.XA.1 {
            self.x_t1d0 = (self.XA.1, self.xt1() * self.t_t1d0(eta, eta1))
        }
        self.x_t1d0.1
    }
    fn x_t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        if self.x_t1d1.0 != self.XA.1 {
            self.x_t1d1 = (
                self.XA.1,
                self.xt2() * self.t_t1d0(eta, eta1) * self.t_t0d1(eta)
                    + self.xt1() * self.t_t1d1(eta, eta1),
            )
        }
        self.x_t1d1.1
    }
    fn x_t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        if self.x_t1d2.0 != self.XA.1 {
            self.x_t1d2 = (
                self.XA.1,
                self.xt3() * self.t_t1d0(eta, eta1) * self.t_t0d1(eta).powi(2)
                    + 2.0 * self.xt2() * self.t_t1d1(eta, eta1) * self.t_t0d1(eta)
                    + self.xt2() * self.t_t1d0(eta, eta1) * self.t_t0d2(eta)
                    + self.xt1() * self.t_t1d2(eta, eta1),
            )
        }
        self.x_t1d2.1
    }
    fn x_t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        if self.x_t1d3.0 != self.XA.1 {
            self.x_t1d3 = (
                self.XA.1,
                self.xt4() * self.t_t1d0(eta, eta1) * self.t_t0d1(eta).powi(3)
                    + 3.0 * self.xt3() * self.t_t1d1(eta, eta1) * self.t_t0d1(eta).powi(2)
                    + 3.0
                        * self.xt3()
                        * self.t_t1d0(eta, eta1)
                        * self.t_t0d1(eta)
                        * self.t_t0d2(eta)
                    + 3.0 * self.xt2() * self.t_t1d2(eta, eta1) * self.t_t0d1(eta)
                    + 3.0 * self.xt2() * self.t_t1d1(eta, eta1) * self.t_t0d2(eta)
                    + self.xt2() * self.t_t1d0(eta, eta1) * self.t_t0d3(eta)
                    + self.xt1() * self.t_t1d3(eta, eta1),
            )
        }
        self.x_t1d3.1
    }
    fn x_t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        if self.x_t2d0.0 != self.XA.1 {
            self.x_t2d0 = (
                self.XA.1,
                self.xt2() * self.t_t1d0(eta, eta1).powi(2)
                    + self.xt1() * self.t_t2d0(eta, eta1, eta2),
            )
        }
        self.x_t2d0.1
    }
}
#[allow(non_snake_case)]
impl PcSaftYglPure {
    fn xt1(&mut self) -> f64 {
        if self.xt1.0 != self.XA.1 {
            let XA = self.XA.1;
            self.xt1 = (
                XA,
                (XA * (XA + self.b_assoc * XA - self.b_assoc)).powi(2)
                    / ((1.0 + self.b_assoc) * XA.powi(2) - 2.0 * (1.0 + self.b_assoc) * XA
                        + self.b_assoc),
            )
        }
        self.xt1.1
    }
    fn xt2(&mut self) -> f64 {
        if self.xt2.0 != self.XA.1 {
            let XA = self.XA.1;
            self.xt2 = (
                XA,
                2.0 * (XA * (XA + self.b_assoc * XA - self.b_assoc)).powi(3)
                    / ((1.0 + self.b_assoc) * XA.powi(2) - 2.0 * (1.0 + self.b_assoc) * XA
                        + self.b_assoc)
                        .powi(3)
                    * ((1.0 + self.b_assoc).powi(2) * XA.powi(3)
                        - 3.0 * (1.0 + self.b_assoc).powi(2) * XA.powi(2)
                        + (3.0 * self.b_assoc * (1.0 + self.b_assoc) * XA - self.b_assoc.powi(2))),
            )
        }
        self.xt2.1
    }
    fn xt3(&mut self) -> f64 {
        if self.xt3.0 != self.XA.1 {
            let XA = self.XA.1;
            self.xt3 = (
                XA,
                6.0 * (XA * (XA + self.b_assoc * XA - self.b_assoc)).powi(4)
                    / ((1.0 + self.b_assoc) * XA.powi(2) - 2.0 * (1.0 + self.b_assoc) * XA
                        + self.b_assoc)
                        .powi(5)
                    * ((1.0 + self.b_assoc).powi(4) * XA.powi(6)
                        - 6.0 * (1.0 + self.b_assoc).powi(4) * XA.powi(5)
                        + 5.0
                            * (3.0 * self.b_assoc + 2.0)
                            * (1.0 + self.b_assoc).powi(3)
                            * XA.powi(4)
                        - 20.0 * self.b_assoc * (1.0 + self.b_assoc).powi(3) * XA.powi(3)
                        + self.b_assoc.powi(2)
                            * (15.0 * self.b_assoc + 16.0)
                            * (1.0 + self.b_assoc)
                            * XA.powi(2)
                        - (6.0 * self.b_assoc.powi(3) * (1.0 + self.b_assoc) * XA
                            - self.b_assoc.powi(4))),
            )
        }
        self.xt3.1
    }
    fn xt4(&mut self) -> f64 {
        if self.xt4.0 != self.XA.1 {
            let XA = self.XA.1;
            self.xt4 = (
                XA,
                24.0 * (XA * (XA + self.b_assoc * XA - self.b_assoc)).powi(5)
                    / ((1.0 + self.b_assoc) * XA.powi(2) - 2.0 * (1.0 + self.b_assoc) * XA
                        + self.b_assoc)
                        .powi(7)
                    * ((1.0 + self.b_assoc).powi(6) * XA.powi(9)
                        - 9.0 * (1.0 + self.b_assoc).powi(6) * XA.powi(8)
                        + (36.0 * self.b_assoc + 29.0) * (1.0 + self.b_assoc).powi(5) * XA.powi(7)
                        - 7.0
                            * (12.0 * self.b_assoc + 5.0)
                            * (1.0 + self.b_assoc).powi(5)
                            * XA.powi(6)
                        + 21.0
                            * self.b_assoc
                            * (6.0 * self.b_assoc + 5.0)
                            * (1.0 + self.b_assoc).powi(4)
                            * XA.powi(5)
                        - self.b_assoc.powi(2)
                            * (126.0 * self.b_assoc.powi(2) + 260.0 * self.b_assoc + 135.0)
                            * (1.0 + self.b_assoc).powi(2)
                            * XA.powi(4)
                        + self.b_assoc.powi(3)
                            * (84.0 * self.b_assoc + 95.0)
                            * (1.0 + self.b_assoc).powi(2)
                            * XA.powi(3)
                        - 3.0
                            * self.b_assoc.powi(4)
                            * (12.0 * self.b_assoc + 13.0)
                            * (1.0 + self.b_assoc)
                            * XA.powi(2)
                        + (9.0 * self.b_assoc.powi(5) * (1.0 + self.b_assoc) * XA
                            - self.b_assoc.powi(6))),
            )
        }
        self.xt4.1
    }
}
impl PcSaftYglPure {
    fn t_t0d0(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d0.0 {
            self.t_t0d0 = (
                eta,
                self.kappa_AB_sigma3_dens * (self.epsilon_AB_temp.exp() - 1.0) * self.gii.t0d0(eta),
            )
        }
        self.t_t0d0.1
    }
    fn t_t0d1(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d1.0 {
            self.t_t0d1 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d1(eta) + self.gii.t0d0(eta)),
            )
        }
        self.t_t0d1.1
    }
    fn t_t0d2(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d2.0 {
            self.t_t0d2 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d2(eta) + 2.0 * self.gii.t0d1(eta)),
            )
        }
        self.t_t0d2.1
    }
    fn t_t0d3(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d3.0 {
            self.t_t0d3 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d3(eta) + 3.0 * self.gii.t0d2(eta)),
            )
        }
        self.t_t0d3.1
    }
    fn t_t0d4(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d4.0 {
            self.t_t0d4 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d4(eta) + 4.0 * self.gii.t0d3(eta)),
            )
        }
        self.t_t0d4.1
    }
    fn t_t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d0.0 {
            self.t_t1d0 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0) * self.gii.t1d0(eta, eta1)
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp) * self.gii.t0d0(eta)),
            )
        }
        self.t_t1d0.1
    }
    fn t_t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d1.0 {
            self.t_t1d1 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii.t1d1(eta, eta1) + self.gii.t1d0(eta, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii.t0d1(eta) + self.gii.t0d0(eta))),
            )
        }
        self.t_t1d1.1
    }
    fn t_t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d2.0 {
            self.t_t1d2 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii.t1d2(eta, eta1) + 2.0 * self.gii.t1d1(eta, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii.t0d2(eta) + 2.0 * self.gii.t0d1(eta))),
            )
        }
        self.t_t1d2.1
    }
    fn t_t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d3.0 {
            self.t_t1d3 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii.t1d3(eta, eta1) + 3.0 * self.gii.t1d2(eta, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii.t0d3(eta) + 3.0 * self.gii.t0d2(eta))),
            )
        }
        self.t_t1d3.1
    }
    fn t_t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        if eta != self.t_t2d0.0 {
            self.t_t2d0 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0) * self.gii.t2d0(eta, eta1, eta2)
                        - 2.0
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * self.gii.t1d0(eta, eta1)
                        + self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * (self.epsilon_AB_temp + 2.0)
                            * self.gii.t0d0(eta)),
            )
        }
        self.t_t2d0.1
    }
}
