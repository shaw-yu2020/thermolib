use super::{CTerm, I1Term, I2Term, PcSaftErr};
use crate::algorithms::{brent_zero, romberg_diff};
use crate::f64consts::{FRAC_NA_1E30, FRAC_PI_2, FRAC_PI_6, PI, R};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::iter::zip;
/// PC-SAFT EOS
/// ```
/// use thermolib::PcSaftPure;
/// let (m, sigma, epsilon) = (2.8611, 2.6826, 205.35); // SO2
/// let mut fluid = PcSaftPure::new_fluid(m, sigma, epsilon);
/// fluid.c_flash().unwrap();
/// assert_eq!(fluid.T().unwrap().round(), 438.0);
/// assert_eq!(fluid.p().unwrap().round(), 9099261.0);
/// assert_eq!(fluid.rho().unwrap().round(), 8079.0);
/// fluid.t_flash(273.15).unwrap();
/// assert_eq!(fluid.p_s().unwrap().round(), 156782.0);
/// assert_eq!(fluid.rho_v().unwrap().round(), 71.0);
/// assert_eq!(fluid.rho_l().unwrap().round(), 22100.0);
/// fluid.tp_flash(273.15, 0.1e6).unwrap();
/// assert_eq!(fluid.rho().unwrap().round(), 45.0);
/// let mut methanol = PcSaftPure::new_fluid(1.5255, 3.23, 188.9);
/// methanol.set_2B_assoc_type(0.035176, 2899.5); // kappa_AB epsilon_AB
/// methanol.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(methanol.rho().unwrap().round(), 24676.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)] // For pyclass hhh
pub struct PcSaftPure {
    m: f64,
    sigma3: f64,
    epsilon_temp: f64,
    epsilon: f64,
    temp: f64,
    rho_num: f64,
    rhov_num: f64,
    rhol_num: f64,
    is_single_phase: bool,
    // changed from temp and rho_num
    eta: f64,
    eta_dt1: f64,
    eta_dt2: f64,
    m2e1s3: f64,
    m2e2s3: f64,
    // association term
    assoc_type: AssocType,
    kappa_AB_sigma3: f64,
    epsilon_AB_temp: f64,
    epsilon_AB: f64,
    XA: f64,
    // modified aly_lee_cv0
    cv_B: f64, // cv_B=B/R-1
    cv_C: f64, // cv_C=C/R
    cv_D: f64, // cv_D=D
    cv_E: f64, // cv_E=E/R
    cv_F: f64, // cv_F=F
    // cached variables
    c: CTerm,   // Cterm
    i1: I1Term, // I1Term
    i2: I2Term, // I2Term
    c_t0d0: (f64, f64),
    c_t0d1: (f64, f64),
    c_t0d2: (f64, f64),
    c_t0d3: (f64, f64),
    c_t0d4: (f64, f64),
    c_t1d0: (f64, f64),
    c_t1d1: (f64, f64),
    c_t1d2: (f64, f64),
    c_t1d3: (f64, f64),
    c_t2d0: (f64, f64),
    c1_t0d0: (f64, f64),
    c1_t0d1: (f64, f64),
    c1_t0d2: (f64, f64),
    c1_t0d3: (f64, f64),
    c1_t0d4: (f64, f64),
    c1_t1d0: (f64, f64),
    c1_t1d1: (f64, f64),
    c1_t1d2: (f64, f64),
    c1_t1d3: (f64, f64),
    c1_t2d0: (f64, f64),
    i2_t0d0: (f64, f64),
    i2_t0d1: (f64, f64),
    i2_t0d2: (f64, f64),
    i2_t0d3: (f64, f64),
    i2_t0d4: (f64, f64),
    i2_t1d0: (f64, f64),
    i2_t1d1: (f64, f64),
    i2_t1d2: (f64, f64),
    i2_t1d3: (f64, f64),
    i2_t2d0: (f64, f64),
    gii_t0d0: (f64, f64),
    gii_t0d1: (f64, f64),
    gii_t0d2: (f64, f64),
    gii_t0d3: (f64, f64),
    gii_t0d4: (f64, f64),
    gii_t1d0: (f64, f64),
    gii_t1d1: (f64, f64),
    gii_t1d2: (f64, f64),
    gii_t1d3: (f64, f64),
    gii_t2d0: (f64, f64),
    x_t0d1: (f64, f64),
    x_t0d2: (f64, f64),
    x_t0d3: (f64, f64),
    x_t0d4: (f64, f64),
    x_t1d0: (f64, f64),
    x_t1d1: (f64, f64),
    x_t1d2: (f64, f64),
    x_t1d3: (f64, f64),
    x_t2d0: (f64, f64),
    xt1: (f64, f64),
    xt2: (f64, f64),
    xt3: (f64, f64),
    xt4: (f64, f64),
}
impl PcSaftPure {
    pub fn new_fluid(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self {
            m,
            epsilon,
            epsilon_temp: 0.0,
            sigma3: sigma.powi(3),
            temp: 1.0,
            rho_num: 1E-10,
            rhov_num: 0.0,
            rhol_num: 0.0,
            is_single_phase: true,
            // changed from temp and rho_num
            eta: 0.0,
            eta_dt1: 0.0,
            eta_dt2: 0.0,
            m2e1s3: 0.0,
            m2e2s3: 0.0,
            // association term
            assoc_type: AssocType::Type0,
            kappa_AB_sigma3: 0.0,
            epsilon_AB_temp: 0.0,
            epsilon_AB: 0.0,
            XA: 1.0,
            // modified aly_lee_cv0
            cv_B: 1.5,
            cv_C: 0.0,
            cv_D: 1.0,
            cv_E: 0.0,
            cv_F: 1.0,
            // cached variables
            c: CTerm::new(m),   // CTerm
            i1: I1Term::new(m), // I1Term
            i2: I2Term::new(m), // I1Term
            c_t0d0: (0.0, 0.0),
            c_t0d1: (0.0, 0.0),
            c_t0d2: (0.0, 0.0),
            c_t0d3: (0.0, 0.0),
            c_t0d4: (0.0, 0.0),
            c_t1d0: (0.0, 0.0),
            c_t1d1: (0.0, 0.0),
            c_t1d2: (0.0, 0.0),
            c_t1d3: (0.0, 0.0),
            c_t2d0: (0.0, 0.0),
            c1_t0d0: (0.0, 0.0),
            c1_t0d1: (0.0, 0.0),
            c1_t0d2: (0.0, 0.0),
            c1_t0d3: (0.0, 0.0),
            c1_t0d4: (0.0, 0.0),
            c1_t1d0: (0.0, 0.0),
            c1_t1d1: (0.0, 0.0),
            c1_t1d2: (0.0, 0.0),
            c1_t1d3: (0.0, 0.0),
            c1_t2d0: (0.0, 0.0),
            i2_t0d0: (0.0, 0.0),
            i2_t0d1: (0.0, 0.0),
            i2_t0d2: (0.0, 0.0),
            i2_t0d3: (0.0, 0.0),
            i2_t0d4: (0.0, 0.0),
            i2_t1d0: (0.0, 0.0),
            i2_t1d1: (0.0, 0.0),
            i2_t1d2: (0.0, 0.0),
            i2_t1d3: (0.0, 0.0),
            i2_t2d0: (0.0, 0.0),
            gii_t0d0: (0.0, 0.0),
            gii_t0d1: (0.0, 0.0),
            gii_t0d2: (0.0, 0.0),
            gii_t0d3: (0.0, 0.0),
            gii_t0d4: (0.0, 0.0),
            gii_t1d0: (0.0, 0.0),
            gii_t1d1: (0.0, 0.0),
            gii_t1d2: (0.0, 0.0),
            gii_t1d3: (0.0, 0.0),
            gii_t2d0: (0.0, 0.0),
            x_t0d1: (0.0, 0.0),
            x_t0d2: (0.0, 0.0),
            x_t0d3: (0.0, 0.0),
            x_t0d4: (0.0, 0.0),
            x_t1d0: (0.0, 0.0),
            x_t1d1: (0.0, 0.0),
            x_t1d2: (0.0, 0.0),
            x_t1d3: (0.0, 0.0),
            x_t2d0: (0.0, 0.0),
            xt1: (0.0, 0.0),
            xt2: (0.0, 0.0),
            xt3: (0.0, 0.0),
            xt4: (0.0, 0.0),
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)] // For pymethods hhh
impl PcSaftPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn py_new(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_1_assoc_type(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type1;
        self.kappa_AB_sigma3 = kappa_AB * self.sigma3;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_2B_assoc_type(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type2B;
        self.kappa_AB_sigma3 = kappa_AB * self.sigma3;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_3B_assoc_type(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type3B;
        self.kappa_AB_sigma3 = kappa_AB * self.sigma3;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_aly_lee_cp0(&mut self, B: f64, C: f64, D: f64, E: f64, F: f64) {
        self.cv_B = B / R - 1.0; // cv_B = B/R -1
        self.cv_C = C / R; // cv_C = C/R
        self.cv_D = D; // cv_D = D
        self.cv_E = E / R; // cv_E = E/R
        self.cv_F = F; // cv_F = F
    }
    pub fn print_derivatives(&mut self) {
        self.check_derivatives(true);
    }
    pub fn td_unchecked(&mut self, temp: f64, dens_mol: f64) {
        self.set_temperature_and_number_density(temp, dens_mol * FRAC_NA_1E30);
        self.is_single_phase = true;
    }
    pub fn c_flash(&mut self) -> anyhow::Result<()> {
        // Iteration from temp_c = 1000 dens_c = 1E-10
        let mut temp_c = 1E3;
        let mut dens_c = 1E-10
            / ((FRAC_PI_6 * self.m * self.sigma3)
                * (1.0 - 0.12 * (-3.0 * self.epsilon / temp_c).exp()).powi(3));
        // Define variables
        let mut p_t0d1 = self.calc_p_t0d1(temp_c, dens_c);
        let mut p_t0d2 = self.calc_p_t0d2(temp_c, dens_c);
        let (mut p_t1d1, mut p_t1d2, mut p_t0d3);
        for _i in 1..100000 {
            p_t1d1 = self.calc_p_t1d1(temp_c, dens_c);
            p_t1d2 = self.calc_p_t1d2(temp_c, dens_c);
            p_t0d3 = self.calc_p_t0d3(temp_c, dens_c);
            temp_c -= (p_t0d1 * p_t0d3 - p_t0d2 * p_t0d2) / (p_t1d1 * p_t0d3 - p_t1d2 * p_t0d2);
            dens_c -= (p_t0d1 * p_t1d2 - p_t0d2 * p_t1d1) / (p_t0d2 * p_t1d2 - p_t0d3 * p_t1d1);
            p_t0d1 = self.calc_p_t0d1(temp_c, dens_c);
            p_t0d2 = self.calc_p_t0d2(temp_c, dens_c);
            if p_t0d1.abs() < 1E3 && p_t0d2.abs() < 1E6 {
                self.set_temperature_and_number_density(temp_c, dens_c);
                self.is_single_phase = true;
                return Ok(());
            }
        }
        Err(anyhow!(PcSaftErr::NotConvForC))
    }
    pub fn t_flash(&mut self, temp: f64) -> anyhow::Result<()> {
        let d3 = self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp()).powi(3);
        // Vapor phase: eta = 1E-10
        let rhov_num_guess = 1E-10 / (FRAC_PI_6 * self.m * d3);
        let mut rhov_num = rhov_num_guess;
        let mut p_t0d1_v = self.calc_p_t0d1(temp, rhov_num);
        for _i in 1..100000 {
            if p_t0d1_v.abs() < 10.0 {
                break;
            } else {
                rhov_num -= p_t0d1_v / self.calc_p_t0d2(temp, rhov_num);
                p_t0d1_v = self.calc_p_t0d1(temp, rhov_num);
            }
        }
        let pv_limit = self.calc_p(temp, rhov_num);
        if pv_limit.is_sign_negative() {
            return Err(anyhow!(PcSaftErr::NotConvForT));
        }
        // Liquid phase: eta = 0.5
        let rhol_num_guess = 0.5 / (FRAC_PI_6 * self.m * d3);
        let mut rhol_num = rhol_num_guess;
        let mut p_t0d1_l = self.calc_p_t0d1(temp, rhol_num);
        for _i in 1..100000 {
            if p_t0d1_l.abs() < 10.0 {
                break;
            } else {
                rhol_num -= p_t0d1_l / self.calc_p_t0d2(temp, rhol_num);
                p_t0d1_l = self.calc_p_t0d1(temp, rhol_num);
            }
        }
        let pl_limit = self.calc_p(temp, rhol_num);
        if pl_limit > pv_limit {
            return Err(anyhow!(PcSaftErr::NotConvForT));
        }
        // Iteration for saturation state
        if rhol_num < rhov_num {
            return Err(anyhow!(PcSaftErr::NotConvForT));
        }
        let rhov_max = rhov_num;
        let lnphi_diff = |p: f64| {
            rhov_num = self.calc_density(temp, p, rhov_num_guess);
            if rhov_num.is_nan() && (p / pv_limit - 1.0).abs() < 1E-2 {
                rhov_num = rhov_max;
            };
            rhol_num = self.calc_density(temp, p, rhol_num_guess);
            self.calc_lnphi(temp, rhov_num) - self.calc_lnphi(temp, rhol_num)
        };
        let ps = brent_zero(lnphi_diff, pv_limit - 1.0, pl_limit.max(1.0));
        if ps.is_nan() {
            Err(anyhow!(PcSaftErr::NotConvForT))
        } else {
            self.is_single_phase = false;
            self.temp = temp;
            self.rhov_num = rhov_num;
            self.rhol_num = rhol_num;
            Ok(())
        }
    }
    pub fn tp_flash(&mut self, temp: f64, pres: f64) -> anyhow::Result<()> {
        let d3 = self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp()).powi(3);
        // Iteration from gas phase: eta = 1E-10
        let rhov_num = self.calc_density(temp, pres, 1E-10 / (FRAC_PI_6 * self.m * d3));
        let lnphi_v = if rhov_num.is_nan() {
            f64::INFINITY
        } else {
            self.calc_lnphi(temp, rhov_num)
        };
        // Iteration from liquid phase: eta = 0.5
        let rhol_num = self.calc_density(temp, pres, 0.5 / (FRAC_PI_6 * self.m * d3));
        let lnphi_l = if rhol_num.is_nan() {
            f64::INFINITY
        } else {
            self.calc_lnphi(temp, rhol_num)
        };
        // Select the correct output
        if lnphi_v.is_infinite() && lnphi_l.is_infinite() {
            Err(anyhow!(PcSaftErr::NotConvForTP))
        } else if lnphi_v.is_infinite() {
            self.set_temperature_and_number_density(temp, rhol_num);
            self.is_single_phase = true;
            Ok(())
        } else if lnphi_l.is_infinite() {
            self.set_temperature_and_number_density(temp, rhov_num);
            self.is_single_phase = true;
            Ok(())
        } else {
            if lnphi_v < lnphi_l {
                self.set_temperature_and_number_density(temp, rhov_num);
            } else {
                self.set_temperature_and_number_density(temp, rhol_num);
            }
            self.is_single_phase = true;
            Ok(())
        }
    }
    pub fn T(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.temp)
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn rho(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.rho_num / FRAC_NA_1E30)
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn p(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_p(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn w(&mut self, molar_mass: f64) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok((self.calc_w2(self.temp, self.rho_num) / molar_mass).sqrt())
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn cv(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_cv(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn cp(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_cp(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn cp_res(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_cp_res(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn h_res(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_h_res(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
        }
    }
    pub fn p_s(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
        } else {
            Ok(self.calc_p(self.temp, self.rhov_num) / 2.0
                + self.calc_p(self.temp, self.rhol_num) / 2.0)
        }
    }
    pub fn T_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
        } else {
            Ok(self.temp)
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
        } else {
            Ok(self.rhov_num / FRAC_NA_1E30)
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
        } else {
            Ok(self.rhol_num / FRAC_NA_1E30)
        }
    }
    pub fn B(&mut self, temp: f64) -> f64 {
        let mut rho_num: f64 = 1E-9;
        let mut val_old: f64 = self.calc_r_t0d1(temp, rho_num) / rho_num;
        let mut val_new: f64 = 0.0;
        loop {
            rho_num /= 10.0;
            if rho_num < 1E-30 {
                break;
            }
            val_new = self.calc_r_t0d1(temp, rho_num) / rho_num;
            if (val_new / val_old - 1.0).abs() < 1E-9 {
                break;
            }
            val_old = val_new;
        }
        val_new * FRAC_NA_1E30
    }
    pub fn C(&mut self, temp: f64) -> f64 {
        let mut rho_num: f64 = 1E-9;
        let mut val_old: f64 = self.calc_r_t0d2(temp, rho_num) / rho_num.powi(2);
        let mut val_new: f64 = 0.0;
        loop {
            rho_num /= 10.0;
            if rho_num < 1E-30 {
                break;
            }
            val_new = self.calc_r_t0d2(temp, rho_num) / rho_num.powi(2);
            if (val_new / val_old - 1.0).abs() < 1E-9 {
                break;
            }
            val_old = val_new;
        }
        val_new * FRAC_NA_1E30.powi(2)
    }
    pub fn D(&mut self, temp: f64) -> f64 {
        let mut rho_num: f64 = 1E-9;
        let mut val_old: f64 = self.calc_r_t0d3(temp, rho_num) / rho_num.powi(3);
        let mut val_new: f64 = 0.0;
        loop {
            rho_num /= 10.0;
            if rho_num < 1E-30 {
                break;
            }
            val_new = self.calc_r_t0d3(temp, rho_num) / rho_num.powi(3);
            if (val_new / val_old - 1.0).abs() < 1E-9 {
                break;
            }
            val_old = val_new;
        }
        val_new * FRAC_NA_1E30.powi(3) / 2.0
    }
    pub fn vec_t_flash_g(&mut self, temp: Vec<f64>) -> Vec<f64> {
        temp.into_iter()
            .map(|t| {
                if self.t_flash(t).is_err() {
                    println!("t_flash_g diverge in {} K", t);
                    f64::NAN
                } else {
                    self.rhov_num
                }
            })
            .collect()
    }
    pub fn vec_t_flash_l(&mut self, temp: Vec<f64>) -> Vec<f64> {
        temp.into_iter()
            .map(|t| {
                if self.t_flash(t).is_err() {
                    println!("t_flash_l diverge in {} K", t);
                    f64::NAN
                } else {
                    self.rhol_num
                }
            })
            .collect()
    }
    pub fn vec_p(&mut self, temp: Vec<f64>, dens_num: Vec<f64>) -> Vec<f64> {
        zip(temp, dens_num)
            .map(|(t, d)| self.calc_p(t, d))
            .collect()
    }
    pub fn vec_tp_flash(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
        zip(temp, pres)
            .map(|(t, p)| {
                if self.tp_flash(t, p).is_err() {
                    println!("tp_flash diverge in {} K {} Pa", t, p);
                    f64::NAN
                } else {
                    self.rho_num
                }
            })
            .collect()
    }
    pub fn vec_tp_flash_g(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
        zip(temp, pres)
            .map(|(t, p)| {
                let d3 = self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp()).powi(3);
                // Iteration from gas phase: eta = 1E-10
                self.calc_density(t, p, 1E-10 / (FRAC_PI_6 * self.m * d3))
            })
            .collect()
    }
    pub fn vec_tp_flash_l(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
        zip(temp, pres)
            .map(|(t, p)| {
                let d3 = self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp()).powi(3);
                // Iteration from gas phase: eta = 0.5
                self.calc_density(t, p, 0.5 / (FRAC_PI_6 * self.m * d3))
            })
            .collect()
    }
    pub fn vec_cp(&mut self, temp: Vec<f64>, dens_num: Vec<f64>) -> Vec<f64> {
        zip(temp, dens_num)
            .map(|(t, d)| self.calc_cp(t, d))
            .collect()
    }
    pub fn vec_w(&mut self, temp: Vec<f64>, dens_num: Vec<f64>, molar_mass: f64) -> Vec<f64> {
        zip(temp, dens_num)
            .map(|(t, d)| (self.calc_w2(t, d) / molar_mass).sqrt()) // m/s
            .collect()
    }
}
const FRAC_RE30_NA: f64 = R / FRAC_NA_1E30;
impl PcSaftPure {
    fn calc_p(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp * rho_num * (1.0 + self.calc_r_t0d1(temp, rho_num))
    }
    fn calc_p_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        (FRAC_RE30_NA * temp)
            * (1.0 + 2.0 * self.calc_r_t0d1(temp, rho_num) + self.calc_r_t0d2(temp, rho_num))
    }
    fn calc_p_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp / rho_num
            * (2.0 * self.calc_r_t0d1(temp, rho_num)
                + 4.0 * self.calc_r_t0d2(temp, rho_num)
                + self.calc_r_t0d3(temp, rho_num))
    }
    fn calc_p_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp / rho_num.powi(2)
            * (6.0 * self.calc_r_t0d2(temp, rho_num)
                + 6.0 * self.calc_r_t0d3(temp, rho_num)
                + self.calc_r_t0d4(temp, rho_num))
    }
    fn calc_p_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        (FRAC_RE30_NA * temp)
            * (1.0
                + 2.0 * self.calc_r_t0d1(temp, rho_num)
                + self.calc_r_t0d2(temp, rho_num)
                + 2.0 * self.calc_r_t1d1(temp, rho_num)
                + self.calc_r_t1d2(temp, rho_num))
    }
    fn calc_p_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp / rho_num
            * (2.0 * self.calc_r_t0d1(temp, rho_num)
                + 4.0 * self.calc_r_t0d2(temp, rho_num)
                + self.calc_r_t0d3(temp, rho_num)
                + 2.0 * self.calc_r_t1d1(temp, rho_num)
                + 4.0 * self.calc_r_t1d2(temp, rho_num)
                + self.calc_r_t1d3(temp, rho_num))
    }
    fn calc_w2(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * temp
            * ((1.0 + 2.0 * self.calc_r_t0d1(temp, rho_num) + self.calc_r_t0d2(temp, rho_num))
                + (1.0 + self.calc_r_t0d1(temp, rho_num) + self.calc_r_t1d1(temp, rho_num)).powi(2)
                    / (self.calc_ideal_cv(temp)
                        - 2.0 * self.calc_r_t1d0(temp, rho_num)
                        - self.calc_r_t2d0(temp, rho_num)))
    }
    fn calc_cv(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * (self.calc_ideal_cv(temp)
            - 2.0 * self.calc_r_t1d0(temp, rho_num)
            - self.calc_r_t2d0(temp, rho_num))
    }
    fn calc_cp(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * (self.calc_ideal_cv(temp)
            - 2.0 * self.calc_r_t1d0(temp, rho_num)
            - self.calc_r_t2d0(temp, rho_num)
            + (1.0 + self.calc_r_t0d1(temp, rho_num) + self.calc_r_t1d1(temp, rho_num)).powi(2)
                / (1.0 + 2.0 * self.calc_r_t0d1(temp, rho_num) + self.calc_r_t0d2(temp, rho_num)))
    }
    fn calc_cp_res(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * ((-1.0 - 2.0 * self.calc_r_t1d0(temp, rho_num) - self.calc_r_t2d0(temp, rho_num))
            + (1.0 + self.calc_r_t0d1(temp, rho_num) + self.calc_r_t1d1(temp, rho_num)).powi(2)
                / (1.0 + 2.0 * self.calc_r_t0d1(temp, rho_num) + self.calc_r_t0d2(temp, rho_num)))
    }
    fn calc_h_res(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * temp * (-self.calc_r_t1d0(temp, rho_num) + self.calc_r_t0d1(temp, rho_num))
    }
    fn calc_lnphi(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.calc_r_t0d0(temp, rho_num) + self.calc_r_t0d1(temp, rho_num)
            - (1.0 + self.calc_r_t0d1(temp, rho_num)).ln()
    }
    fn calc_density(&mut self, temp: f64, p: f64, rho_num_guess: f64) -> f64 {
        let mut rho_num = rho_num_guess;
        let (mut p_diff, mut val_p_t0d1, mut rho_num_diff);
        for _i in 1..10000 {
            p_diff = self.calc_p(temp, rho_num) - p;
            if p_diff.abs() < f64::EPSILON {
                return rho_num;
            }
            val_p_t0d1 = self.calc_p_t0d1(temp, rho_num);
            if val_p_t0d1.is_sign_negative() {
                return f64::NAN;
            }
            rho_num_diff = p_diff / val_p_t0d1;
            if rho_num_diff.abs() < f64::EPSILON {
                return rho_num;
            }
            rho_num -= rho_num_diff;
            if rho_num.is_sign_negative() {
                return f64::NAN;
            }
        }
        f64::NAN
    }
}
impl PcSaftPure {
    fn set_temperature_and_number_density(&mut self, temp: f64, rho_num: f64) {
        if temp != self.temp {
            self.temp = temp;
            self.epsilon_temp = self.epsilon / temp;
            self.epsilon_AB_temp = self.epsilon_AB / temp;
            self.m2e1s3 = self.m.powi(2) * self.epsilon_temp * self.sigma3;
            self.m2e2s3 = self.m.powi(2) * self.epsilon_temp.powi(2) * self.sigma3;
        } else if rho_num != self.rho_num {
        } else {
            return;
        }
        self.rho_num = rho_num;
        let d = 1.0 - 0.12 * (-3.0 * self.epsilon_temp).exp();
        let d1 = -0.36 * (-3.0 * self.epsilon_temp).exp() * self.epsilon / temp;
        let d2 = d1 * (3.0 * self.epsilon_temp - 2.0);
        self.eta = FRAC_PI_6 * rho_num * self.m * self.sigma3 * d.powi(3);
        self.eta_dt1 = FRAC_PI_2 * rho_num * self.m * self.sigma3 * d.powi(2) * d1;
        self.eta_dt2 = PI * rho_num * self.m * self.sigma3 * d * d1.powi(2)
            + FRAC_PI_2 * rho_num * self.m * self.sigma3 * d.powi(2) * d2;
        match self.assoc_type {
            AssocType::Type0 => (),
            AssocType::Type1 | AssocType::Type2B => {
                let t = self.rho_num * self.kappa_AB_sigma3 * self.t_t0d0();
                self.XA = (-1.0 + (1.0 + 4.0 * t).sqrt()) / (2.0 * t);
            }
            AssocType::Type3B => {
                let t = self.rho_num * self.kappa_AB_sigma3 * self.t_t0d0();
                self.XA = (-(1.0 - t) + ((1.0 - t).powi(2) + 8.0 * t).sqrt()) / (4.0 * t);
            }
        }
    }
    fn calc_ideal_cv(&mut self, temp: f64) -> f64 {
        self.cv_B
            + self.cv_C * (self.cv_D / temp / (self.cv_D / temp).sinh()).powi(2)
            + self.cv_E * (self.cv_F / temp / (self.cv_F / temp).cosh()).powi(2)
    }
    fn calc_r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d0() + self.disp_t0d0() + self.assoc_t0d0()
    }
    fn calc_r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d1() + self.disp_t0d1() + self.assoc_t0d1()
    }
    fn calc_r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d2() + self.disp_t0d2() + self.assoc_t0d2()
    }
    fn calc_r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d3() + self.disp_t0d3() + self.assoc_t0d3()
    }
    fn calc_r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d4() + self.disp_t0d4() + self.assoc_t0d4()
    }
    fn calc_r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d0() + self.disp_t1d0() + self.assoc_t1d0()
    }
    fn calc_r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d1() + self.disp_t1d1() + self.assoc_t1d1()
    }
    fn calc_r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d2() + self.disp_t1d2() + self.assoc_t1d2()
    }
    fn calc_r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d3() + self.disp_t1d3() + self.assoc_t1d3()
    }
    fn calc_r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t2d0() + self.disp_t2d0() + self.assoc_t2d0()
    }
}
impl PcSaftPure {
    fn hc_t0d0(&mut self) -> f64 {
        self.m * self.hs_t0d0() + (1.0 - self.m) * self.gii_t0d0().ln()
    }
    fn hc_t0d1(&mut self) -> f64 {
        self.m * self.hs_t0d1() + (1.0 - self.m) / self.gii_t0d0() * self.gii_t0d1()
    }
    fn hc_t0d2(&mut self) -> f64 {
        self.m * self.hs_t0d2()
            + (1.0 - self.m)
                * (-1.0 / self.gii_t0d0().powi(2) * self.gii_t0d1().powi(2)
                    + 1.0 / self.gii_t0d0() * self.gii_t0d2())
    }
    fn hc_t0d3(&mut self) -> f64 {
        self.m * self.hs_t0d3()
            + (1.0 - self.m)
                * (2.0 * self.gii_t0d1().powi(3) / self.gii_t0d0().powi(3)
                    - 3.0 * self.gii_t0d1() * self.gii_t0d2() / self.gii_t0d0().powi(2)
                    + self.gii_t0d3() / self.gii_t0d0())
    }
    fn hc_t0d4(&mut self) -> f64 {
        self.m * self.hs_t0d4()
            + (1.0 - self.m)
                * (-6.0 / self.gii_t0d0().powi(4) * self.gii_t0d1().powi(4)
                    + 12.0 / self.gii_t0d0().powi(3) * self.gii_t0d1().powi(2) * self.gii_t0d2()
                    - 4.0 / self.gii_t0d0().powi(2) * self.gii_t0d1() * self.gii_t0d3()
                    - 3.0 / self.gii_t0d0().powi(2) * self.gii_t0d2().powi(2)
                    + self.gii_t0d4() / self.gii_t0d0())
    }
    fn hc_t1d0(&mut self) -> f64 {
        self.m * self.hs_t1d0() + (1.0 - self.m) / self.gii_t0d0() * self.gii_t1d0()
    }
    fn hc_t1d1(&mut self) -> f64 {
        self.m * self.hs_t1d1()
            + (1.0 - self.m)
                * (-self.gii_t1d0() * self.gii_t0d1() / self.gii_t0d0().powi(2)
                    + self.gii_t1d1() / self.gii_t0d0())
    }
    fn hc_t1d2(&mut self) -> f64 {
        self.m * self.hs_t1d2()
            + (1.0 - self.m)
                * (2.0 / self.gii_t0d0().powi(3) * self.gii_t1d0() * self.gii_t0d1().powi(2)
                    - 2.0 / self.gii_t0d0().powi(2) * self.gii_t1d1() * self.gii_t0d1()
                    - self.gii_t1d0() * self.gii_t0d2() / self.gii_t0d0().powi(2)
                    + self.gii_t1d2() / self.gii_t0d0())
    }
    fn hc_t1d3(&mut self) -> f64 {
        self.m * self.hs_t1d3()
            + (1.0 - self.m)
                * (-6.0 / self.gii_t0d0().powi(4) * self.gii_t1d0() * self.gii_t0d1().powi(3)
                    + 6.0 / self.gii_t0d0().powi(3) * self.gii_t1d1() * self.gii_t0d1().powi(2)
                    + 6.0 / self.gii_t0d0().powi(3)
                        * self.gii_t1d0()
                        * self.gii_t0d1()
                        * self.gii_t0d2()
                    - 3.0 / self.gii_t0d0().powi(2) * self.gii_t1d1() * self.gii_t0d2()
                    - 3.0 / self.gii_t0d0().powi(2) * self.gii_t1d2() * self.gii_t0d1()
                    - self.gii_t1d0() * self.gii_t0d3() / self.gii_t0d0().powi(2)
                    + self.gii_t1d3() / self.gii_t0d0())
    }
    fn hc_t2d0(&mut self) -> f64 {
        self.m * self.hs_t2d0()
            + (1.0 - self.m)
                * (-self.gii_t1d0().powi(2) / self.gii_t0d0().powi(2)
                    + self.gii_t2d0() / self.gii_t0d0())
    }
}
impl PcSaftPure {
    fn hs_t0d0(&self) -> f64 {
        self.eta * (4.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(2)
    }
    fn hs_t0d1(&self) -> f64 {
        self.eta * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
    }
    fn hs_t0d2(&self) -> f64 {
        self.eta.powi(2) * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
    fn hs_t0d3(&self) -> f64 {
        self.eta.powi(3) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn hs_t0d4(&self) -> f64 {
        self.eta.powi(4) * (168.0 - 48.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn hs_t1d0(&self) -> f64 {
        self.eta_dt1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
    }
    fn hs_t1d1(&self) -> f64 {
        self.eta_dt1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta_dt1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
    fn hs_t1d2(&self) -> f64 {
        2.0 * self.eta_dt1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
            + self.eta_dt1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn hs_t1d3(&self) -> f64 {
        3.0 * self.eta_dt1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
            + self.eta_dt1 * self.eta.powi(3) * (168.0 - 48.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn hs_t2d0(&self) -> f64 {
        self.eta_dt2 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta_dt1.powi(2) * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
}
impl PcSaftPure {
    fn gii_t0d0(&mut self) -> f64 {
        if self.gii_t0d0.0 != self.eta {
            self.gii_t0d0 = (self.eta, (1.0 - 0.5 * self.eta) / (1.0 - self.eta).powi(3))
        }
        self.gii_t0d0.1
    }
    fn gii_t0d1(&mut self) -> f64 {
        if self.gii_t0d1.0 != self.eta {
            self.gii_t0d1 = (
                self.eta,
                self.eta * (2.5 - self.eta) / (1.0 - self.eta).powi(4),
            )
        }
        self.gii_t0d1.1
    }
    fn gii_t0d2(&mut self) -> f64 {
        if self.gii_t0d2.0 != self.eta {
            self.gii_t0d2 = (
                self.eta,
                self.eta.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5),
            )
        }
        self.gii_t0d2.1
    }
    fn gii_t0d3(&mut self) -> f64 {
        if self.gii_t0d3.0 != self.eta {
            self.gii_t0d3 = (
                self.eta,
                self.eta.powi(3) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6),
            )
        }
        self.gii_t0d3.1
    }
    fn gii_t0d4(&mut self) -> f64 {
        if self.gii_t0d4.0 != self.eta {
            self.gii_t0d4 = (
                self.eta,
                self.eta.powi(4) * (240.0 - 60.0 * self.eta) / (1.0 - self.eta).powi(7),
            )
        }
        self.gii_t0d4.1
    }
    fn gii_t1d0(&mut self) -> f64 {
        if self.gii_t1d0.0 != self.eta {
            self.gii_t1d0 = (
                self.eta,
                self.eta_dt1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4),
            )
        }
        self.gii_t1d0.1
    }
    fn gii_t1d1(&mut self) -> f64 {
        if self.gii_t1d1.0 != self.eta {
            self.gii_t1d1 = (
                self.eta,
                self.eta_dt1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
                    + self.eta_dt1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5),
            )
        }
        self.gii_t1d1.1
    }
    fn gii_t1d2(&mut self) -> f64 {
        if self.gii_t1d2.0 != self.eta {
            self.gii_t1d2 = (
                self.eta,
                2.0 * self.eta_dt1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
                    + self.eta_dt1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta)
                        / (1.0 - self.eta).powi(6),
            )
        }
        self.gii_t1d2.1
    }
    fn gii_t1d3(&mut self) -> f64 {
        if self.gii_t1d3.0 != self.eta {
            self.gii_t1d3 = (
                self.eta,
                3.0 * self.eta_dt1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta)
                    / (1.0 - self.eta).powi(6)
                    + self.eta_dt1 * self.eta.powi(3) * (240.0 - 60.0 * self.eta)
                        / (1.0 - self.eta).powi(7),
            )
        }
        self.gii_t1d3.1
    }
    fn gii_t2d0(&mut self) -> f64 {
        if self.gii_t2d0.0 != self.eta {
            self.gii_t2d0 = (
                self.eta,
                self.eta_dt2 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
                    + self.eta_dt1.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5),
            )
        }
        self.gii_t2d0.1
    }
}
impl PcSaftPure {
    fn disp_t0d0(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * self.i1_t0d0(self.eta)
                + self.m * self.m2e2s3 * self.c1_t0d0() * self.i2_t0d0())
    }
    fn disp_t0d1(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (self.i1_t0d0(self.eta) + self.i1_t0d1(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1_t0d0() * self.i2_t0d0()
                        + self.c1_t0d1() * self.i2_t0d0()
                        + self.c1_t0d0() * self.i2_t0d1()))
    }
    fn disp_t0d2(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (2.0 * self.i1_t0d1(self.eta) + self.i1_t0d2(self.eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * self.c1_t0d1() * self.i2_t0d0()
                        + 2.0 * self.c1_t0d0() * self.i2_t0d1()
                        + self.c1_t0d2() * self.i2_t0d0()
                        + 2.0 * self.c1_t0d1() * self.i2_t0d1()
                        + self.c1_t0d0() * self.i2_t0d2()))
    }
    fn disp_t0d3(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (3.0 * self.i1_t0d2(self.eta) + self.i1_t0d3(self.eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * self.c1_t0d2() * self.i2_t0d0()
                        + 6.0 * self.c1_t0d1() * self.i2_t0d1()
                        + 3.0 * self.c1_t0d0() * self.i2_t0d2()
                        + self.c1_t0d3() * self.i2_t0d0()
                        + 3.0 * self.c1_t0d2() * self.i2_t0d1()
                        + 3.0 * self.c1_t0d1() * self.i2_t0d2()
                        + self.c1_t0d0() * self.i2_t0d3()))
    }
    fn disp_t0d4(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (4.0 * self.i1_t0d3(self.eta) + self.i1_t0d4(self.eta))
                + self.m
                    * self.m2e2s3
                    * (4.0 * self.c1_t0d3() * self.i2_t0d0()
                        + 12.0 * self.c1_t0d2() * self.i2_t0d1()
                        + 12.0 * self.c1_t0d1() * self.i2_t0d2()
                        + 4.0 * self.c1_t0d0() * self.i2_t0d3()
                        + self.c1_t0d4() * self.i2_t0d0()
                        + 4.0 * self.c1_t0d3() * self.i2_t0d1()
                        + 6.0 * self.c1_t0d2() * self.i2_t0d2()
                        + 4.0 * self.c1_t0d1() * self.i2_t0d3()
                        + self.c1_t0d0() * self.i2_t0d4()))
    }
    fn disp_t1d0(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (self.i1_t1d0(self.eta) - self.i1_t0d0(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1_t1d0() * self.i2_t0d0() + self.c1_t0d0() * self.i2_t1d0()
                        - 2.0 * self.c1_t0d0() * self.i2_t0d0()))
    }
    fn disp_t1d1(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (self.i1_t1d0(self.eta) - self.i1_t0d0(self.eta) + self.i1_t1d1(self.eta)
                    - self.i1_t0d1(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1_t1d0() * self.i2_t0d0() + self.c1_t0d0() * self.i2_t1d0()
                        - 2.0 * self.c1_t0d0() * self.i2_t0d0()
                        + self.c1_t1d1() * self.i2_t0d0()
                        + self.c1_t0d1() * self.i2_t1d0()
                        - 2.0 * self.c1_t0d1() * self.i2_t0d0()
                        + self.c1_t1d0() * self.i2_t0d1()
                        + self.c1_t0d0() * self.i2_t1d1()
                        - 2.0 * self.c1_t0d0() * self.i2_t0d1()))
    }
    fn disp_t1d2(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (2.0 * self.i1_t1d1(self.eta) - 2.0 * self.i1_t0d1(self.eta)
                    + self.i1_t1d2(self.eta)
                    - self.i1_t0d2(self.eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * self.c1_t1d1() * self.i2_t0d0()
                        + 2.0 * self.c1_t0d1() * self.i2_t1d0()
                        - 4.0 * self.c1_t0d1() * self.i2_t0d0()
                        + 2.0 * self.c1_t1d0() * self.i2_t0d1()
                        + 2.0 * self.c1_t0d0() * self.i2_t1d1()
                        - 4.0 * self.c1_t0d0() * self.i2_t0d1()
                        + self.c1_t1d2() * self.i2_t0d0()
                        + self.c1_t0d2() * self.i2_t1d0()
                        - 2.0 * self.c1_t0d2() * self.i2_t0d0()
                        + 2.0 * self.c1_t1d1() * self.i2_t0d1()
                        + 2.0 * self.c1_t0d1() * self.i2_t1d1()
                        - 4.0 * self.c1_t0d1() * self.i2_t0d1()
                        + self.c1_t1d0() * self.i2_t0d2()
                        + self.c1_t0d0() * self.i2_t1d2()
                        - 2.0 * self.c1_t0d0() * self.i2_t0d2()))
    }
    fn disp_t1d3(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (3.0 * self.i1_t1d2(self.eta) - 3.0 * self.i1_t0d2(self.eta)
                    + self.i1_t1d3(self.eta)
                    - self.i1_t0d3(self.eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * self.c1_t1d2() * self.i2_t0d0()
                        + 3.0 * self.c1_t0d2() * self.i2_t1d0()
                        - 6.0 * self.c1_t0d2() * self.i2_t0d0()
                        + 6.0 * self.c1_t1d1() * self.i2_t0d1()
                        + 6.0 * self.c1_t0d1() * self.i2_t1d1()
                        - 12.0 * self.c1_t0d1() * self.i2_t0d1()
                        + 3.0 * self.c1_t1d0() * self.i2_t0d2()
                        + 3.0 * self.c1_t0d0() * self.i2_t1d2()
                        - 6.0 * self.c1_t0d0() * self.i2_t0d2()
                        + self.c1_t1d3() * self.i2_t0d0()
                        + self.c1_t0d3() * self.i2_t1d0()
                        - 2.0 * self.c1_t0d3() * self.i2_t0d0()
                        + 3.0 * self.c1_t1d2() * self.i2_t0d1()
                        + 3.0 * self.c1_t0d2() * self.i2_t1d1()
                        - 6.0 * self.c1_t0d2() * self.i2_t0d1()
                        + 3.0 * self.c1_t1d1() * self.i2_t0d2()
                        + 3.0 * self.c1_t0d1() * self.i2_t1d2()
                        - 6.0 * self.c1_t0d1() * self.i2_t0d2()
                        + self.c1_t1d0() * self.i2_t0d3()
                        + self.c1_t0d0() * self.i2_t1d3()
                        - 2.0 * self.c1_t0d0() * self.i2_t0d3()))
    }
    fn disp_t2d0(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (self.i1_t2d0(self.eta) + 2.0 * self.i1_t0d0(self.eta)
                    - 2.0 * self.i1_t1d0(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1_t2d0() * self.i2_t0d0()
                        + self.c1_t0d0() * self.i2_t2d0()
                        + 6.0 * self.c1_t0d0() * self.i2_t0d0()
                        + 2.0 * self.c1_t1d0() * self.i2_t1d0()
                        - 4.0 * self.c1_t1d0() * self.i2_t0d0()
                        - 4.0 * self.c1_t0d0() * self.i2_t1d0()))
    }
}
impl PcSaftPure {
    fn c1_t0d0(&mut self) -> f64 {
        if self.c1_t0d0.0 != self.eta {
            self.c1_t0d0 = (self.eta, self.c_t0d0().recip())
        }
        self.c1_t0d0.1
    }
    fn c1_t0d1(&mut self) -> f64 {
        if self.c1_t0d1.0 != self.eta {
            self.c1_t0d1 = (self.eta, -self.c_t0d1() / self.c_t0d0().powi(2))
        }
        self.c1_t0d1.1
    }
    fn c1_t0d2(&mut self) -> f64 {
        if self.c1_t0d2.0 != self.eta {
            self.c1_t0d2 = (
                self.eta,
                2.0 * self.c_t0d1().powi(2) / self.c_t0d0().powi(3)
                    - self.c_t0d2() / self.c_t0d0().powi(2),
            )
        }
        self.c1_t0d2.1
    }
    fn c1_t0d3(&mut self) -> f64 {
        if self.c1_t0d3.0 != self.eta {
            self.c1_t0d3 = (
                self.eta,
                -6.0 * self.c_t0d1().powi(3) / self.c_t0d0().powi(4)
                    + 6.0 * self.c_t0d1() * self.c_t0d2() / self.c_t0d0().powi(3)
                    - self.c_t0d3() / self.c_t0d0().powi(2),
            )
        }
        self.c1_t0d3.1
    }
    fn c1_t0d4(&mut self) -> f64 {
        if self.c1_t0d4.0 != self.eta {
            self.c1_t0d4 = (
                self.eta,
                24.0 * self.c_t0d1().powi(4) / self.c_t0d0().powi(5)
                    - 36.0 * self.c_t0d1().powi(2) * self.c_t0d2() / self.c_t0d0().powi(4)
                    + 8.0 * self.c_t0d1() * self.c_t0d3() / self.c_t0d0().powi(3)
                    + 6.0 * self.c_t0d2().powi(2) / self.c_t0d0().powi(3)
                    - self.c_t0d4() / self.c_t0d0().powi(2),
            )
        }
        self.c1_t0d4.1
    }
    fn c1_t1d0(&mut self) -> f64 {
        if self.c1_t1d0.0 != self.eta {
            self.c1_t1d0 = (self.eta, -self.c_t1d0() / self.c_t0d0().powi(2))
        }
        self.c1_t1d0.1
    }
    fn c1_t1d1(&mut self) -> f64 {
        if self.c1_t1d1.0 != self.eta {
            self.c1_t1d1 = (
                self.eta,
                2.0 * self.c_t1d0() * self.c_t0d1() / self.c_t0d0().powi(3)
                    - self.c_t1d1() / self.c_t0d0().powi(2),
            )
        }
        self.c1_t1d1.1
    }
    fn c1_t1d2(&mut self) -> f64 {
        if self.c1_t1d2.0 != self.eta {
            self.c1_t1d2 = (
                self.eta,
                -6.0 * self.c_t1d0() * self.c_t0d1().powi(2) / self.c_t0d0().powi(4)
                    + 4.0 * self.c_t1d1() * self.c_t0d1() / self.c_t0d0().powi(3)
                    + 2.0 * self.c_t1d0() * self.c_t0d2() / self.c_t0d0().powi(3)
                    - self.c_t1d2() / self.c_t0d0().powi(2),
            )
        }
        self.c1_t1d2.1
    }
    fn c1_t1d3(&mut self) -> f64 {
        if self.c1_t1d3.0 != self.eta {
            self.c1_t1d3 = (
                self.eta,
                24.0 * self.c_t1d0() * self.c_t0d1().powi(3) / self.c_t0d0().powi(5)
                    - 18.0 * self.c_t1d1() * self.c_t0d1().powi(2) / self.c_t0d0().powi(4)
                    - 18.0 * self.c_t1d0() * self.c_t0d1() * self.c_t0d2() / self.c_t0d0().powi(4)
                    + 6.0 * self.c_t1d1() * self.c_t0d2() / self.c_t0d0().powi(3)
                    + 6.0 * self.c_t1d2() * self.c_t0d1() / self.c_t0d0().powi(3)
                    + 2.0 * self.c_t1d0() * self.c_t0d3() / self.c_t0d0().powi(3)
                    - self.c_t1d3() / self.c_t0d0().powi(2),
            )
        }
        self.c1_t1d3.1
    }
    fn c1_t2d0(&mut self) -> f64 {
        if self.c1_t2d0.0 != self.eta {
            self.c1_t2d0 = (
                self.eta,
                2.0 * self.c_t1d0().powi(2) / self.c_t0d0().powi(3)
                    - self.c_t2d0() / self.c_t0d0().powi(2),
            )
        }
        self.c1_t2d0.1
    }
}
impl PcSaftPure {
    fn c_t0d0(&mut self) -> f64 {
        if self.c_t0d0.0 != self.eta {
            self.c_t0d0 = (self.eta, self.c.eta0(self.eta))
        }
        self.c_t0d0.1
    }
    fn c_t0d1(&mut self) -> f64 {
        if self.c_t0d1.0 != self.eta {
            self.c_t0d1 = (self.eta, self.eta * self.c.eta1(self.eta))
        }
        self.c_t0d1.1
    }
    fn c_t0d2(&mut self) -> f64 {
        if self.c_t0d2.0 != self.eta {
            self.c_t0d2 = (self.eta, self.eta.powi(2) * self.c.eta2(self.eta))
        }
        self.c_t0d2.1
    }
    fn c_t0d3(&mut self) -> f64 {
        if self.c_t0d3.0 != self.eta {
            self.c_t0d3 = (self.eta, self.eta.powi(3) * self.c.eta3(self.eta))
        }
        self.c_t0d3.1
    }
    fn c_t0d4(&mut self) -> f64 {
        if self.c_t0d4.0 != self.eta {
            self.c_t0d4 = (self.eta, self.eta.powi(4) * self.c.eta4(self.eta))
        }
        self.c_t0d4.1
    }
    fn c_t1d0(&mut self) -> f64 {
        if self.c_t1d0.0 != self.eta {
            self.c_t1d0 = (self.eta, self.eta_dt1 * self.c.eta1(self.eta))
        }
        self.c_t1d0.1
    }
    fn c_t1d1(&mut self) -> f64 {
        if self.c_t1d1.0 != self.eta {
            self.c_t1d1 = (
                self.eta,
                self.eta_dt1 * (self.c.eta1(self.eta) + self.eta * self.c.eta2(self.eta)),
            )
        }
        self.c_t1d1.1
    }
    fn c_t1d2(&mut self) -> f64 {
        if self.c_t1d2.0 != self.eta {
            self.c_t1d2 = (
                self.eta,
                (self.eta_dt1 * self.eta)
                    * (2.0 * self.c.eta2(self.eta) + self.eta * self.c.eta3(self.eta)),
            )
        }
        self.c_t1d2.1
    }
    fn c_t1d3(&mut self) -> f64 {
        if self.c_t1d3.0 != self.eta {
            self.c_t1d3 = (
                self.eta,
                (self.eta_dt1 * self.eta.powi(2))
                    * (3.0 * self.c.eta3(self.eta) + self.eta * self.c.eta4(self.eta)),
            )
        }
        self.c_t1d3.1
    }
    fn c_t2d0(&mut self) -> f64 {
        if self.c_t2d0.0 != self.eta {
            self.c_t2d0 = (
                self.eta,
                self.eta_dt2 * self.c.eta1(self.eta) + self.eta_dt1.powi(2) * self.c.eta2(self.eta),
            )
        }
        self.c_t2d0.1
    }
}
impl PcSaftPure {
    fn i1_t0d0(&mut self, eta: f64) -> f64 {
        self.i1.eta0(eta)
    }
    fn i1_t0d1(&mut self, eta: f64) -> f64 {
        eta * self.i1.eta1(eta)
    }
    fn i1_t0d2(&mut self, eta: f64) -> f64 {
        eta.powi(2) * self.i1.eta2(eta)
    }
    fn i1_t0d3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.i1.eta3(eta)
    }
    fn i1_t0d4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.i1.eta4(eta)
    }
    fn i1_t1d0(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.i1.eta1(eta)
    }
    fn i1_t1d1(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * (self.i1.eta1(eta) + self.eta * self.i1.eta2(eta))
    }
    fn i1_t1d2(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta * (2.0 * self.i1.eta2(eta) + self.eta * self.i1.eta3(eta))
    }
    fn i1_t1d3(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta.powi(2) * (3.0 * self.i1.eta3(eta) + self.eta * self.i1.eta4(eta))
    }
    fn i1_t2d0(&mut self, eta: f64) -> f64 {
        self.eta_dt2 * self.i1.eta1(eta) + self.eta_dt1.powi(2) * self.i1.eta2(eta)
    }
}
impl PcSaftPure {
    fn i2_t0d0(&mut self) -> f64 {
        if self.i2_t0d0.0 != self.eta {
            self.i2_t0d0 = (self.eta, self.i2.eta0(self.eta))
        }
        self.i2_t0d0.1
    }
    fn i2_t0d1(&mut self) -> f64 {
        if self.i2_t0d1.0 != self.eta {
            self.i2_t0d1 = (self.eta, self.eta * self.i2.eta1(self.eta))
        }
        self.i2_t0d1.1
    }
    fn i2_t0d2(&mut self) -> f64 {
        if self.i2_t0d2.0 != self.eta {
            self.i2_t0d2 = (self.eta, self.eta.powi(2) * self.i2.eta2(self.eta))
        }
        self.i2_t0d2.1
    }
    fn i2_t0d3(&mut self) -> f64 {
        if self.i2_t0d3.0 != self.eta {
            self.i2_t0d3 = (self.eta, self.eta.powi(3) * self.i2.eta3(self.eta))
        }
        self.i2_t0d3.1
    }
    fn i2_t0d4(&mut self) -> f64 {
        if self.i2_t0d4.0 != self.eta {
            self.i2_t0d4 = (self.eta, self.eta.powi(4) * self.i2.eta4(self.eta))
        }
        self.i2_t0d4.1
    }
    fn i2_t1d0(&mut self) -> f64 {
        if self.i2_t1d0.0 != self.eta {
            self.i2_t1d0 = (self.eta, self.eta_dt1 * self.i2.eta1(self.eta))
        }
        self.i2_t1d0.1
    }
    fn i2_t1d1(&mut self) -> f64 {
        if self.i2_t1d1.0 != self.eta {
            self.i2_t1d1 = (
                self.eta,
                self.eta_dt1 * (self.i2.eta1(self.eta) + self.eta * self.i2.eta2(self.eta)),
            )
        }
        self.i2_t1d1.1
    }
    fn i2_t1d2(&mut self) -> f64 {
        if self.i2_t1d2.0 != self.eta {
            self.i2_t1d2 = (
                self.eta,
                (self.eta_dt1 * self.eta)
                    * (2.0 * self.i2.eta2(self.eta) + self.eta * self.i2.eta3(self.eta)),
            )
        }
        self.i2_t1d2.1
    }
    fn i2_t1d3(&mut self) -> f64 {
        if self.i2_t1d3.0 != self.eta {
            self.i2_t1d3 = (
                self.eta,
                (self.eta_dt1 * self.eta.powi(2))
                    * (3.0 * self.i2.eta3(self.eta) + self.eta * self.i2.eta4(self.eta)),
            )
        }
        self.i2_t1d3.1
    }
    fn i2_t2d0(&mut self) -> f64 {
        if self.i2_t2d0.0 != self.eta {
            self.i2_t2d0 = (
                self.eta,
                self.eta_dt2 * self.i2.eta1(self.eta)
                    + self.eta_dt1.powi(2) * self.i2.eta2(self.eta),
            )
        }
        self.i2_t2d0.1
    }
}
enum AssocType {
    Type0,
    Type1,
    Type2B,
    Type3B,
}
impl PcSaftPure {
    fn assoc_t0d0(&self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.XA.ln() - self.XA / 2.0 + 0.5,
            AssocType::Type2B => 2.0 * self.XA.ln() - self.XA + 1.0,
            AssocType::Type3B => {
                2.0 * self.XA.ln() + (2.0 * self.XA - 1.0).ln() - 2.0 * self.XA + 2.0
            }
        }
    }
    fn assoc_t0d1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t0d1::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t0d1::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d1::<1>(self.XA) + self.site_t0d1::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t0d2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t0d2::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t0d2::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d2::<1>(self.XA) + self.site_t0d2::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t0d3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t0d3::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t0d3::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d3::<1>(self.XA) + self.site_t0d3::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t0d4(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t0d4::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t0d4::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d4::<1>(self.XA) + self.site_t0d4::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t1d0::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t1d0::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d0::<1>(self.XA) + self.site_t1d0::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t1d1::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t1d1::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d1::<1>(self.XA) + self.site_t1d1::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t1d2::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t1d2::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d2::<1>(self.XA) + self.site_t1d2::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t1d3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t1d3::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t1d3::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t1d3::<1>(self.XA) + self.site_t1d3::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
    fn assoc_t2d0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.site_t2d0::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t2d0::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t2d0::<1>(self.XA) + self.site_t2d0::<2>(2.0 * self.XA - 1.0)
            }
        }
    }
}
impl PcSaftPure {
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
impl PcSaftPure {
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
impl PcSaftPure {
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
impl PcSaftPure {
    fn xt1(&mut self) -> f64 {
        if self.xt1.0 != self.XA {
            self.xt1 = (
                self.XA,
                (self.rho_num * self.kappa_AB_sigma3)
                    * match self.assoc_type {
                        AssocType::Type1 | AssocType::Type2B => self.XA.powi(3) / (self.XA - 2.0),
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
                2.0 * (self.rho_num * self.kappa_AB_sigma3).powi(2)
                    * match self.assoc_type {
                        AssocType::Type1 | AssocType::Type2B => {
                            self.XA.powi(5) / (self.XA - 2.0).powi(3) * (self.XA - 3.0)
                        }
                        AssocType::Type3B => {
                            (self.XA * (2.0 * self.XA - 1.0)).powi(3)
                                / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(3)
                                * (4.0 * self.XA.powi(3) - 12.0 * self.XA.powi(2) + 6.0 * self.XA
                                    - 1.0)
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
                6.0 * (self.rho_num * self.kappa_AB_sigma3).powi(3)
                    * match self.assoc_type {
                        AssocType::Type1 | AssocType::Type2B => {
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
                24.0 * (self.rho_num * self.kappa_AB_sigma3).powi(4)
                    * match self.assoc_type {
                        AssocType::Type1 | AssocType::Type2B => {
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
                        _ => 0.0,
                    },
            )
        }
        self.xt4.1
    }
}
impl PcSaftPure {
    fn t_t0d0(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * self.gii_t0d0()
    }
    fn t_t0d1(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * (self.gii_t0d1() + self.gii_t0d0())
    }
    fn t_t0d2(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * (self.gii_t0d2() + 2.0 * self.gii_t0d1())
    }
    fn t_t0d3(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * (self.gii_t0d3() + 3.0 * self.gii_t0d2())
    }
    fn t_t0d4(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * (self.gii_t0d4() + 4.0 * self.gii_t0d3())
    }
    fn t_t1d0(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * self.gii_t1d0()
            - self.epsilon_AB_temp.exp() * self.epsilon_AB_temp * self.gii_t0d0()
    }
    fn t_t1d1(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * (self.gii_t1d1() + self.gii_t1d0())
            - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                * (self.gii_t0d1() + self.gii_t0d0())
    }
    fn t_t1d2(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * (self.gii_t1d2() + 2.0 * self.gii_t1d1())
            - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                * (self.gii_t0d2() + 2.0 * self.gii_t0d1())
    }
    fn t_t1d3(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * (self.gii_t1d3() + 3.0 * self.gii_t1d2())
            - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                * (self.gii_t0d3() + 3.0 * self.gii_t0d2())
    }
    fn t_t2d0(&mut self) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0) * self.gii_t2d0()
            - 2.0 * self.epsilon_AB_temp.exp() * self.epsilon_AB_temp * self.gii_t1d0()
            + (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                * (self.epsilon_AB_temp + 2.0)
                * self.gii_t0d0()
    }
}
impl PcSaftPure {
    pub fn check_derivatives(&mut self, print_val: bool) {
        let (t, d) = (self.temp, self.rho_num);
        if print_val {
            println!("[rT0D0 == rT0D0] calc_r_t0d0 ={}", self.calc_r_t0d0(t, d));
        }
        let compare_val = |val_calc: f64, val_diff: f64| {
            assert_eq!(
                &val_calc.abs().to_string()[0..11],
                &val_diff.abs().to_string()[0..11]
            )
        };
        // derivative for density
        let val_calc = self.calc_r_t0d1(t, d) / d;
        let val_diff = romberg_diff(|dx: f64| self.calc_r_t0d0(t, dx), d);
        if print_val {
            println!("[rT0D1 == rT0D1] calc_r_t0d1 ={}", val_calc);
            println!("[rT0D0 -> rT0D1] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density
        let val_calc = self.calc_r_t0d2(t, d) / d.powi(2);
        let val_diff = romberg_diff(|dx: f64| self.calc_r_t0d1(t, dx) / dx, d);
        if print_val {
            println!("[rT0D2 == rT0D2] calc_r_t0d2 ={}", val_calc);
            println!("[rT0D1 -> rT0D2] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density
        let val_calc = self.calc_r_t0d3(t, d) / d.powi(3);
        let val_diff = romberg_diff(|dx: f64| self.calc_r_t0d2(t, dx) / dx.powi(2), d);
        if print_val {
            println!("[rT0D3 == rT0D3] calc_r_t0d3 ={}", val_calc);
            println!("[rT0D2 -> rT0D3] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density+density
        let val_calc = self.calc_r_t0d4(t, d) / d.powi(4);
        let val_diff = romberg_diff(|dx: f64| self.calc_r_t0d3(t, dx) / dx.powi(3), d);
        if print_val {
            println!("[rT0D4 == rT0D4] calc_r_t0d4 ={}", val_calc);
            println!("[rT0D3 -> rT0D4] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature
        let val_calc = self.calc_r_t1d0(t, d) / t;
        let val_diff = romberg_diff(|tx: f64| self.calc_r_t0d0(tx, d), t);
        if print_val {
            println!("[rT1D0 == rT1D0] calc_r_t1d0 ={}", val_calc);
            println!("[rT0D0 -> rT1D0] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density
        let val_calc = self.calc_r_t1d1(t, d) / t / d;
        let val_diff = romberg_diff(|dx: f64| self.calc_r_t1d0(t, dx) / t, d);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_r_t1d1 ={}", val_calc);
            println!("[rT1D0 -> rT1D1] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.calc_r_t0d1(tx, d) / d, t);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_r_t1d1 ={}", val_calc);
            println!("[rT0D1 -> rT1D1] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density
        let val_calc = self.calc_r_t1d2(t, d) / t / d.powi(2);
        let val_diff = romberg_diff(|dx: f64| self.calc_r_t1d1(t, dx) / t / dx, d);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_r_t1d2 ={}", val_calc);
            println!("[rT1D1 -> rT1D2] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.calc_r_t0d2(tx, d) / d.powi(2), t);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_r_t1d2 ={}", val_calc);
            println!("[rT0D2 -> rT1D2] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density+density
        let val_calc = self.calc_r_t1d3(t, d) / t / d.powi(3);
        let val_diff = romberg_diff(|dx: f64| self.calc_r_t1d2(t, dx) / t / dx.powi(2), d);
        if print_val {
            println!("[rT1D3 == rT1D3] calc_r_t1d3 ={}", val_calc);
            println!("[rT1D2 -> rT1D3] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.calc_r_t0d3(tx, d) / d.powi(3), t);
        if print_val {
            println!("[rT1D3 == rT1D3] calc_r_t1d3 ={}", val_calc);
            println!("[rT0D3 -> rT1D3] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+temperature
        let val_calc = self.calc_r_t2d0(t, d) / t.powi(2);
        let val_diff = romberg_diff(|tx: f64| self.calc_r_t1d0(tx, d) / tx, t);
        if print_val {
            println!("[rT2D0 == rT2D0] calc_r_t2d0 ={}", val_calc);
            println!("[rT1D0 -> rT2D0] romberg_dif ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
    }
}
/// unit test
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pc_saft_pure() {
        let (m, sigma, epsilon) = (2.8611, 2.6826, 205.35); // SO2
        let mut fluid = PcSaftPure::new_fluid(m, sigma, epsilon);
        fluid.c_flash().unwrap();
        fluid.check_derivatives(false);
        let temp_min = (0.6 * fluid.T().unwrap()).floor() as u32;
        let temp_max = fluid.T().unwrap().ceil() as u32;
        let (mut p_s, mut rho_v, mut rho_l): (f64, f64, f64) = (0.0, 0.0, f64::INFINITY);
        for temp in temp_min..temp_max {
            fluid.t_flash(temp as f64).unwrap();
            if fluid.p_s().unwrap() < p_s
                || fluid.rho_v().unwrap() < rho_v
                || fluid.rho_l().unwrap() > rho_l
            {
                panic!()
            } else {
                p_s = fluid.p_s().unwrap();
                rho_v = fluid.rho_v().unwrap();
                rho_l = fluid.rho_l().unwrap();
            }
        }
        let mut methanol = PcSaftPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_2B_assoc_type(0.035176, 2899.5);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
        let mut methanol = PcSaftPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_3B_assoc_type(0.035176, 2899.5);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
    }
}
