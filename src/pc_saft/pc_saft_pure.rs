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

/*
以下是旧代码
*/

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
                let t = self.tT0D0();
                self.XA = (-1.0 + (1.0 + 4.0 * t).sqrt()) / (2.0 * t);
            }
            AssocType::Type3B => {
                let t = self.tT0D0();
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
        self.hcT0D0() + self.dispT0D0() + self.assocT0D0()
    }
    fn calc_r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D1() + self.dispT0D1() + self.assocT0D1()
    }
    fn calc_r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D2() + self.dispT0D2() + self.assocT0D2()
    }
    fn calc_r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D3() + self.dispT0D3() + self.assocT0D3()
    }
    fn calc_r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D4() + self.dispT0D4() + self.assocT0D4()
    }
    fn calc_r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D0() + self.dispT1D0() + self.assocT1D0()
    }
    fn calc_r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D1() + self.dispT1D1() + self.assocT1D1()
    }
    fn calc_r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D2() + self.dispT1D2() + self.assocT1D2()
    }
    fn calc_r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D3() + self.dispT1D3() + self.assocT1D3()
    }
    fn calc_r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT2D0() + self.dispT2D0() + self.assocT2D0()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn hcT0D0(&mut self) -> f64 {
        self.m * self.hsT0D0() + (1.0 - self.m) * self.giiT0D0().ln()
    }
    fn hcT0D1(&mut self) -> f64 {
        self.m * self.hsT0D1() + (1.0 - self.m) / self.giiT0D0() * self.giiT0D1()
    }
    fn hcT0D2(&mut self) -> f64 {
        self.m * self.hsT0D2()
            + (1.0 - self.m)
                * (-1.0 / self.giiT0D0().powi(2) * self.giiT0D1().powi(2)
                    + 1.0 / self.giiT0D0() * self.giiT0D2())
    }
    fn hcT0D3(&mut self) -> f64 {
        self.m * self.hsT0D3()
            + (1.0 - self.m)
                * (2.0 * self.giiT0D1().powi(3) / self.giiT0D0().powi(3)
                    - 3.0 * self.giiT0D1() * self.giiT0D2() / self.giiT0D0().powi(2)
                    + self.giiT0D3() / self.giiT0D0())
    }
    fn hcT0D4(&mut self) -> f64 {
        self.m * self.hsT0D4()
            + (1.0 - self.m)
                * (-6.0 / self.giiT0D0().powi(4) * self.giiT0D1().powi(4)
                    + 12.0 / self.giiT0D0().powi(3) * self.giiT0D1().powi(2) * self.giiT0D2()
                    - 4.0 / self.giiT0D0().powi(2) * self.giiT0D1() * self.giiT0D3()
                    - 3.0 / self.giiT0D0().powi(2) * self.giiT0D2().powi(2)
                    + self.giiT0D4() / self.giiT0D0())
    }
    fn hcT1D0(&mut self) -> f64 {
        self.m * self.hsT1D0() + (1.0 - self.m) / self.giiT0D0() * self.giiT1D0()
    }
    fn hcT1D1(&mut self) -> f64 {
        self.m * self.hsT1D1()
            + (1.0 - self.m)
                * (-self.giiT1D0() * self.giiT0D1() / self.giiT0D0().powi(2)
                    + self.giiT1D1() / self.giiT0D0())
    }
    fn hcT1D2(&mut self) -> f64 {
        let giiT1D0: f64 = self.giiT1D0();
        self.m * self.hsT1D2()
            + (1.0 - self.m)
                * (2.0 / self.giiT0D0().powi(3) * giiT1D0 * self.giiT0D1().powi(2)
                    - 2.0 / self.giiT0D0().powi(2) * self.giiT1D1() * self.giiT0D1()
                    - giiT1D0 * self.giiT0D2() / self.giiT0D0().powi(2)
                    + self.giiT1D2() / self.giiT0D0())
    }
    fn hcT1D3(&mut self) -> f64 {
        let giiT1D0: f64 = self.giiT1D0();
        let giiT1D1: f64 = self.giiT1D1();
        self.m * self.hsT1D3()
            + (1.0 - self.m)
                * (-6.0 / self.giiT0D0().powi(4) * giiT1D0 * self.giiT0D1().powi(3)
                    + 6.0 / self.giiT0D0().powi(3) * giiT1D1 * self.giiT0D1().powi(2)
                    + 6.0 / self.giiT0D0().powi(3) * giiT1D0 * self.giiT0D1() * self.giiT0D2()
                    - 3.0 / self.giiT0D0().powi(2) * giiT1D1 * self.giiT0D2()
                    - 3.0 / self.giiT0D0().powi(2) * self.giiT1D2() * self.giiT0D1()
                    - giiT1D0 * self.giiT0D3() / self.giiT0D0().powi(2)
                    + self.giiT1D3() / self.giiT0D0())
    }
    fn hcT2D0(&mut self) -> f64 {
        self.m * self.hsT2D0()
            + (1.0 - self.m)
                * (-self.giiT1D0().powi(2) / self.giiT0D0().powi(2)
                    + self.giiT2D0() / self.giiT0D0())
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn hsT0D0(&self) -> f64 {
        self.eta * (4.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(2)
    }
    fn hsT0D1(&self) -> f64 {
        self.eta * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
    }
    fn hsT0D2(&self) -> f64 {
        self.eta.powi(2) * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
    fn hsT0D3(&self) -> f64 {
        self.eta.powi(3) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn hsT0D4(&self) -> f64 {
        self.eta.powi(4) * (168.0 - 48.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn hsT1D0(&self) -> f64 {
        self.eta_dt1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
    }
    fn hsT1D1(&self) -> f64 {
        self.eta_dt1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta_dt1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
    fn hsT1D2(&self) -> f64 {
        2.0 * self.eta_dt1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
            + self.eta_dt1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn hsT1D3(&self) -> f64 {
        3.0 * self.eta_dt1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
            + self.eta_dt1 * self.eta.powi(3) * (168.0 - 48.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn hsT2D0(&self) -> f64 {
        self.eta_dt2 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta_dt1.powi(2) * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn giiT0D0(&mut self) -> f64 {
        (1.0 - 0.5 * self.eta) / (1.0 - self.eta).powi(3)
    }
    fn giiT0D1(&mut self) -> f64 {
        self.eta * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
    }
    fn giiT0D2(&mut self) -> f64 {
        self.eta.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn giiT0D3(&self) -> f64 {
        self.eta.powi(3) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn giiT0D4(&self) -> f64 {
        self.eta.powi(4) * (240.0 - 60.0 * self.eta) / (1.0 - self.eta).powi(7)
    }
    fn giiT1D0(&self) -> f64 {
        self.eta_dt1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
    }
    fn giiT1D1(&self) -> f64 {
        self.eta_dt1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
            + self.eta_dt1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn giiT1D2(&self) -> f64 {
        2.0 * self.eta_dt1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
            + self.eta_dt1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn giiT1D3(&self) -> f64 {
        3.0 * self.eta_dt1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
            + self.eta_dt1 * self.eta.powi(3) * (240.0 - 60.0 * self.eta) / (1.0 - self.eta).powi(7)
    }
    fn giiT2D0(&self) -> f64 {
        self.eta_dt2 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
            + self.eta_dt1.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn dispT0D0(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * self.i1T0D0(self.eta)
                + self.m * self.m2e2s3 * self.c1T0D0() * self.i2T0D0())
    }
    fn dispT0D1(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (self.i1T0D0(self.eta) + self.i1T0D1(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1T0D0() * self.i2T0D0()
                        + self.c1T0D1() * self.i2T0D0()
                        + self.c1T0D0() * self.i2T0D1()))
    }
    fn dispT0D2(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (2.0 * self.i1T0D1(self.eta) + self.i1T0D2(self.eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * self.c1T0D1() * self.i2T0D0()
                        + 2.0 * self.c1T0D0() * self.i2T0D1()
                        + self.c1T0D2() * self.i2T0D0()
                        + 2.0 * self.c1T0D1() * self.i2T0D1()
                        + self.c1T0D0() * self.i2T0D2()))
    }
    fn dispT0D3(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (3.0 * self.i1T0D2(self.eta) + self.i1T0D3(self.eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * self.c1T0D2() * self.i2T0D0()
                        + 6.0 * self.c1T0D1() * self.i2T0D1()
                        + 3.0 * self.c1T0D0() * self.i2T0D2()
                        + self.c1T0D3(self.eta) * self.i2T0D0()
                        + 3.0 * self.c1T0D2() * self.i2T0D1()
                        + 3.0 * self.c1T0D1() * self.i2T0D2()
                        + self.c1T0D0() * self.i2T0D3(self.eta)))
    }
    fn dispT0D4(&mut self) -> f64 {
        let c1T0D3 = self.c1T0D3(self.eta);
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (4.0 * self.i1T0D3(self.eta) + self.i1T0D4(self.eta))
                + self.m
                    * self.m2e2s3
                    * (4.0 * c1T0D3 * self.i2T0D0()
                        + 12.0 * self.c1T0D2() * self.i2T0D1()
                        + 12.0 * self.c1T0D1() * self.i2T0D2()
                        + 4.0 * self.c1T0D0() * self.i2T0D3(self.eta)
                        + self.c1T0D4(self.eta) * self.i2T0D0()
                        + 4.0 * c1T0D3 * self.i2T0D1()
                        + 6.0 * self.c1T0D2() * self.i2T0D2()
                        + 4.0 * self.c1T0D1() * self.i2T0D3(self.eta)
                        + self.c1T0D0() * self.i2T0D4(self.eta)))
    }
    fn dispT1D0(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (self.i1T1D0(self.eta) - self.i1T0D0(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1T1D0(self.eta) * self.i2T0D0()
                        + self.c1T0D0() * self.i2T1D0(self.eta)
                        - 2.0 * self.c1T0D0() * self.i2T0D0()))
    }
    fn dispT1D1(&mut self) -> f64 {
        let c1T1D0 = self.c1T1D0(self.eta);
        let i2T1D0 = self.i2T1D0(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (self.i1T1D0(self.eta) - self.i1T0D0(self.eta) + self.i1T1D1(self.eta)
                    - self.i1T0D1(self.eta))
                + self.m
                    * self.m2e2s3
                    * (c1T1D0 * self.i2T0D0() + self.c1T0D0() * i2T1D0
                        - 2.0 * self.c1T0D0() * self.i2T0D0()
                        + self.c1T1D1(self.eta) * self.i2T0D0()
                        + self.c1T0D1() * i2T1D0
                        - 2.0 * self.c1T0D1() * self.i2T0D0()
                        + c1T1D0 * self.i2T0D1()
                        + self.c1T0D0() * self.i2T1D1(self.eta)
                        - 2.0 * self.c1T0D0() * self.i2T0D1()))
    }
    fn dispT1D2(&mut self) -> f64 {
        let c1T1D0 = self.c1T1D0(self.eta);
        let c1T1D1 = self.c1T1D1(self.eta);
        let i2T1D0 = self.i2T1D0(self.eta);
        let i2T1D1 = self.i2T1D1(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (2.0 * self.i1T1D1(self.eta) - 2.0 * self.i1T0D1(self.eta)
                    + self.i1T1D2(self.eta)
                    - self.i1T0D2(self.eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * c1T1D1 * self.i2T0D0() + 2.0 * self.c1T0D1() * i2T1D0
                        - 4.0 * self.c1T0D1() * self.i2T0D0()
                        + 2.0 * c1T1D0 * self.i2T0D1()
                        + 2.0 * self.c1T0D0() * i2T1D1
                        - 4.0 * self.c1T0D0() * self.i2T0D1()
                        + self.c1T1D2(self.eta) * self.i2T0D0()
                        + self.c1T0D2() * i2T1D0
                        - 2.0 * self.c1T0D2() * self.i2T0D0()
                        + 2.0 * c1T1D1 * self.i2T0D1()
                        + 2.0 * self.c1T0D1() * i2T1D1
                        - 4.0 * self.c1T0D1() * self.i2T0D1()
                        + c1T1D0 * self.i2T0D2()
                        + self.c1T0D0() * self.i2T1D2(self.eta)
                        - 2.0 * self.c1T0D0() * self.i2T0D2()))
    }
    fn dispT1D3(&mut self) -> f64 {
        let c1T0D3 = self.c1T0D3(self.eta);
        let c1T1D0 = self.c1T1D0(self.eta);
        let c1T1D1 = self.c1T1D1(self.eta);
        let c1T1D2 = self.c1T1D2(self.eta);
        let i2T1D0 = self.i2T1D0(self.eta);
        let i2T1D1 = self.i2T1D1(self.eta);
        let i2T1D2 = self.i2T1D2(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (3.0 * self.i1T1D2(self.eta) - 3.0 * self.i1T0D2(self.eta)
                    + self.i1T1D3(self.eta)
                    - self.i1T0D3(self.eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * self.c1T1D2(self.eta) * self.i2T0D0() + 3.0 * self.c1T0D2() * i2T1D0
                        - 6.0 * self.c1T0D2() * self.i2T0D0()
                        + 6.0 * c1T1D1 * self.i2T0D1()
                        + 6.0 * self.c1T0D1() * i2T1D1
                        - 12.0 * self.c1T0D1() * self.i2T0D1()
                        + 3.0 * c1T1D0 * self.i2T0D2()
                        + 3.0 * self.c1T0D0() * self.i2T1D2(self.eta)
                        - 6.0 * self.c1T0D0() * self.i2T0D2()
                        + self.c1T1D3(self.eta) * self.i2T0D0()
                        + c1T0D3 * i2T1D0
                        - 2.0 * c1T0D3 * self.i2T0D0()
                        + 3.0 * c1T1D2 * self.i2T0D1()
                        + 3.0 * self.c1T0D2() * i2T1D1
                        - 6.0 * self.c1T0D2() * self.i2T0D1()
                        + 3.0 * c1T1D1 * self.i2T0D2()
                        + 3.0 * self.c1T0D1() * i2T1D2
                        - 6.0 * self.c1T0D1() * self.i2T0D2()
                        + c1T1D0 * self.i2T0D3(self.eta)
                        + self.c1T0D0() * self.i2T1D3(self.eta)
                        - 2.0 * self.c1T0D0() * self.i2T0D3(self.eta)))
    }
    fn dispT2D0(&mut self) -> f64 {
        let c1T1D0 = self.c1T1D0(self.eta);
        let i2T1D0 = self.i2T1D0(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (self.i1T2D0(self.eta) + 2.0 * self.i1T0D0(self.eta)
                    - 2.0 * self.i1T1D0(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1T2D0(self.eta) * self.i2T0D0()
                        + self.c1T0D0() * self.i2T2D0(self.eta)
                        + 6.0 * self.c1T0D0() * self.i2T0D0()
                        + 2.0 * c1T1D0 * i2T1D0
                        - 4.0 * c1T1D0 * self.i2T0D0()
                        - 4.0 * self.c1T0D0() * i2T1D0))
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn c1T0D0(&mut self) -> f64 {
        self.cT0D0().recip()
    }
    fn c1T0D1(&mut self) -> f64 {
        -self.cT0D1() / self.cT0D0().powi(2)
    }
    fn c1T0D2(&mut self) -> f64 {
        2.0 * self.cT0D1().powi(2) / self.cT0D0().powi(3) - self.cT0D2() / self.cT0D0().powi(2)
    }
    fn c1T0D3(&mut self, eta: f64) -> f64 {
        -6.0 * self.cT0D1().powi(3) / self.cT0D0().powi(4)
            + 6.0 * self.cT0D1() * self.cT0D2() / self.cT0D0().powi(3)
            - self.cT0D3(eta) / self.cT0D0().powi(2)
    }
    fn c1T0D4(&mut self, eta: f64) -> f64 {
        24.0 * self.cT0D1().powi(4) / self.cT0D0().powi(5)
            - 36.0 * self.cT0D1().powi(2) * self.cT0D2() / self.cT0D0().powi(4)
            + 8.0 * self.cT0D1() * self.cT0D3(eta) / self.cT0D0().powi(3)
            + 6.0 * self.cT0D2().powi(2) / self.cT0D0().powi(3)
            - self.cT0D4(eta) / self.cT0D0().powi(2)
    }
    fn c1T1D0(&mut self, eta: f64) -> f64 {
        -self.cT1D0(eta) / self.cT0D0().powi(2)
    }
    fn c1T1D1(&mut self, eta: f64) -> f64 {
        2.0 * self.cT1D0(eta) * self.cT0D1() / self.cT0D0().powi(3)
            - self.cT1D1(eta) / self.cT0D0().powi(2)
    }
    fn c1T1D2(&mut self, eta: f64) -> f64 {
        let cT1D0 = self.cT1D0(eta);
        -6.0 * cT1D0 * self.cT0D1().powi(2) / self.cT0D0().powi(4)
            + 4.0 * self.cT1D1(eta) * self.cT0D1() / self.cT0D0().powi(3)
            + 2.0 * cT1D0 * self.cT0D2() / self.cT0D0().powi(3)
            - self.cT1D2(eta) / self.cT0D0().powi(2)
    }
    fn c1T1D3(&mut self, eta: f64) -> f64 {
        let cT1D0 = self.cT1D0(eta);
        let cT1D1 = self.cT1D1(eta);
        24.0 * cT1D0 * self.cT0D1().powi(3) / self.cT0D0().powi(5)
            - 18.0 * cT1D1 * self.cT0D1().powi(2) / self.cT0D0().powi(4)
            - 18.0 * cT1D0 * self.cT0D1() * self.cT0D2() / self.cT0D0().powi(4)
            + 6.0 * cT1D1 * self.cT0D2() / self.cT0D0().powi(3)
            + 6.0 * self.cT1D2(eta) * self.cT0D1() / self.cT0D0().powi(3)
            + 2.0 * cT1D0 * self.cT0D3(eta) / self.cT0D0().powi(3)
            - self.cT1D3(eta) / self.cT0D0().powi(2)
    }
    fn c1T2D0(&mut self, eta: f64) -> f64 {
        2.0 * self.cT1D0(eta).powi(2) / self.cT0D0().powi(3)
            - self.cT2D0(eta) / self.cT0D0().powi(2)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn cT0D0(&mut self) -> f64 {
        self.c.eta0(self.eta)
    }
    fn cT0D1(&mut self) -> f64 {
        self.eta * self.c.eta1(self.eta)
    }
    fn cT0D2(&mut self) -> f64 {
        self.eta.powi(2) * self.c.eta2(self.eta)
    }
    fn cT0D3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.c.eta3(eta)
    }
    fn cT0D4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.c.eta4(eta)
    }
    fn cT1D0(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.c.eta1(eta)
    }
    fn cT1D1(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * (self.c.eta1(eta) + self.eta * self.c.eta2(eta))
    }
    fn cT1D2(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta * (2.0 * self.c.eta2(eta) + self.eta * self.c.eta3(eta))
    }
    fn cT1D3(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta.powi(2) * (3.0 * self.c.eta3(eta) + self.eta * self.c.eta4(eta))
    }
    fn cT2D0(&mut self, eta: f64) -> f64 {
        self.eta_dt2 * self.c.eta1(eta) + self.eta_dt1.powi(2) * self.c.eta2(eta)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn i1T0D0(&mut self, eta: f64) -> f64 {
        self.i1.eta0(eta)
    }
    fn i1T0D1(&mut self, eta: f64) -> f64 {
        eta * self.i1.eta1(eta)
    }
    fn i1T0D2(&mut self, eta: f64) -> f64 {
        eta.powi(2) * self.i1.eta2(eta)
    }
    fn i1T0D3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.i1.eta3(eta)
    }
    fn i1T0D4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.i1.eta4(eta)
    }
    fn i1T1D0(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.i1.eta1(eta)
    }
    fn i1T1D1(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * (self.i1.eta1(eta) + self.eta * self.i1.eta2(eta))
    }
    fn i1T1D2(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta * (2.0 * self.i1.eta2(eta) + self.eta * self.i1.eta3(eta))
    }
    fn i1T1D3(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta.powi(2) * (3.0 * self.i1.eta3(eta) + self.eta * self.i1.eta4(eta))
    }
    fn i1T2D0(&mut self, eta: f64) -> f64 {
        self.eta_dt2 * self.i1.eta1(eta) + self.eta_dt1.powi(2) * self.i1.eta2(eta)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn i2T0D0(&mut self) -> f64 {
        self.i2.eta0(self.eta)
    }
    fn i2T0D1(&mut self) -> f64 {
        self.eta * self.i2.eta1(self.eta)
    }
    fn i2T0D2(&mut self) -> f64 {
        self.eta.powi(2) * self.i2.eta2(self.eta)
    }
    fn i2T0D3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.i2.eta3(eta)
    }
    fn i2T0D4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.i2.eta4(eta)
    }
    fn i2T1D0(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.i2.eta1(eta)
    }
    fn i2T1D1(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * (self.i2.eta1(eta) + self.eta * self.i2.eta2(eta))
    }
    fn i2T1D2(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta * (2.0 * self.i2.eta2(eta) + self.eta * self.i2.eta3(eta))
    }
    fn i2T1D3(&mut self, eta: f64) -> f64 {
        self.eta_dt1 * self.eta.powi(2) * (3.0 * self.i2.eta3(eta) + self.eta * self.i2.eta4(eta))
    }
    fn i2T2D0(&mut self, eta: f64) -> f64 {
        self.eta_dt2 * self.i2.eta1(eta) + self.eta_dt1.powi(2) * self.i2.eta2(eta)
    }
}
enum AssocType {
    Type0,
    Type1,
    Type2B,
    Type3B,
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn assocT0D0(&self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.XA.ln() - self.XA / 2.0 + 0.5,
            AssocType::Type2B => 2.0 * self.XA.ln() - self.XA + 1.0,
            AssocType::Type3B => {
                2.0 * self.XA.ln() + (2.0 * self.XA - 1.0).ln() - 2.0 * self.XA + 2.0
            }
        }
    }
    fn assocT0D1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT0D1(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT0D1(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT0D1(1.0, self.XA) + self.siteT0D1(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT0D2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT0D2(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT0D2(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT0D2(1.0, self.XA) + self.siteT0D2(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT0D3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT0D3(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT0D3(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT0D3(1.0, self.XA) + self.siteT0D3(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT0D4(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT0D4(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT0D4(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT0D4(1.0, self.XA) + self.siteT0D4(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT1D0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT1D0(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT1D0(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT1D0(1.0, self.XA) + self.siteT1D0(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT1D1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT1D1(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT1D1(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT1D1(1.0, self.XA) + self.siteT1D1(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT1D2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT1D2(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT1D2(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT1D2(1.0, self.XA) + self.siteT1D2(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT1D3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT1D3(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT1D3(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT1D3(1.0, self.XA) + self.siteT1D3(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
    fn assocT2D0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 => self.siteT2D0(1.0, self.XA),
            AssocType::Type2B => 2.0 * self.siteT2D0(1.0, self.XA),
            AssocType::Type3B => {
                2.0 * self.siteT2D0(1.0, self.XA) + self.siteT2D0(2.0, 2.0 * self.XA - 1.0)
            }
        }
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn siteT0D1(&mut self, c: f64, X: f64) -> f64 {
        c * self.siteX1(X) * self.XT0D1()
    }
    fn siteT0D2(&mut self, c: f64, X: f64) -> f64 {
        c.powi(2) * self.siteX2(X) * self.XT0D1().powi(2) + c * self.siteX1(X) * self.XT0D2()
    }
    fn siteT0D3(&mut self, c: f64, X: f64) -> f64 {
        c.powi(3) * self.siteX3(X) * self.XT0D1().powi(3)
            + 3.0 * c.powi(2) * self.siteX2(X) * self.XT0D1() * self.XT0D2()
            + c * self.siteX1(X) * self.XT0D3()
    }
    fn siteT0D4(&mut self, c: f64, X: f64) -> f64 {
        c.powi(4) * self.siteX4(X) * self.XT0D1().powi(4)
            + 6.0 * c.powi(3) * self.siteX3(X) * self.XT0D1().powi(2) * self.XT0D2()
            + 3.0 * c.powi(2) * self.siteX2(X) * self.XT0D2().powi(2)
            + 4.0 * c.powi(2) * self.siteX2(X) * self.XT0D1() * self.XT0D3()
            + c * self.siteX1(X) * self.XT0D4()
    }
    fn siteT1D0(&mut self, c: f64, X: f64) -> f64 {
        c * self.siteX1(X) * self.XT1D0()
    }
    fn siteT1D1(&mut self, c: f64, X: f64) -> f64 {
        c.powi(2) * self.siteX2(X) * self.XT1D0() * self.XT0D1() + c * self.siteX1(X) * self.XT1D1()
    }
    fn siteT1D2(&mut self, c: f64, X: f64) -> f64 {
        c.powi(3) * self.siteX3(X) * self.XT1D0() * self.XT0D1().powi(2)
            + 2.0 * c.powi(2) * self.siteX2(X) * self.XT1D1() * self.XT0D1()
            + c.powi(2) * self.siteX2(X) * self.XT1D0() * self.XT0D2()
            + c * self.siteX1(X) * self.XT1D2()
    }
    fn siteT1D3(&mut self, c: f64, X: f64) -> f64 {
        c.powi(4) * self.siteX4(X) * self.XT1D0() * self.XT0D1().powi(3)
            + 3.0 * c.powi(3) * self.siteX3(X) * self.XT1D1() * self.XT0D1().powi(2)
            + 3.0 * c.powi(3) * self.siteX3(X) * self.XT1D0() * self.XT0D1() * self.XT0D2()
            + 3.0 * c.powi(2) * self.siteX2(X) * self.XT1D2() * self.XT0D1()
            + 3.0 * c.powi(2) * self.siteX2(X) * self.XT1D1() * self.XT0D2()
            + c.powi(2) * self.siteX2(X) * self.XT1D0() * self.XT0D3()
            + c * self.siteX1(X) * self.XT1D3()
    }
    fn siteT2D0(&mut self, c: f64, X: f64) -> f64 {
        c.powi(2) * self.siteX2(X) * self.XT1D0().powi(2) + c * self.siteX1(X) * self.XT2D0()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn siteX1(&self, X: f64) -> f64 {
        1.0 / X - 1.0 / 2.0
    }
    fn siteX2(&self, X: f64) -> f64 {
        -1.0 / X.powi(2)
    }
    fn siteX3(&self, X: f64) -> f64 {
        2.0 / X.powi(3)
    }
    fn siteX4(&self, X: f64) -> f64 {
        -6.0 / X.powi(4)
    }
}
impl PcSaftPure {
    fn xt1(&self) -> f64 {
        match self.assoc_type {
            AssocType::Type1 | AssocType::Type2B => self.XA.powi(3) / (self.XA - 2.0),
            AssocType::Type3B => {
                (self.XA * (2.0 * self.XA - 1.0)).powi(2)
                    / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0)
            }
            AssocType::Type0 => 0.0,
        }
    }
    fn xt2(&self) -> f64 {
        2.0 * match self.assoc_type {
            AssocType::Type1 | AssocType::Type2B => {
                self.XA.powi(5) / (self.XA - 2.0).powi(3) * (self.XA - 3.0)
            }
            AssocType::Type3B => {
                (self.XA * (2.0 * self.XA - 1.0)).powi(3)
                    / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(3)
                    * (4.0 * self.XA.powi(3) - 12.0 * self.XA.powi(2) + 6.0 * self.XA - 1.0)
            }
            AssocType::Type0 => 0.0,
        }
    }
    fn xt3(&self) -> f64 {
        6.0 * match self.assoc_type {
            AssocType::Type1 | AssocType::Type2B => {
                self.XA.powi(7) / (self.XA - 2.0).powi(5) * (self.XA.powi(2) - 6.0 * self.XA + 10.0)
            }
            AssocType::Type3B => {
                (self.XA * (2.0 * self.XA - 1.0)).powi(4)
                    / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(5)
                    * (16.0 * self.XA.powi(6) - 96.0 * self.XA.powi(5)
                        + (200.0 * self.XA.powi(4) - 160.0 * self.XA.powi(3))
                        + (62.0 * self.XA.powi(2) - 12.0 * self.XA + 1.0))
            }
            AssocType::Type0 => 0.0,
        }
    }
    fn xt4(&self) -> f64 {
        24.0 * match self.assoc_type {
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
        }
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn XT0D1(&mut self) -> f64 {
        self.xt1() * self.tT0D1()
    }
    fn XT0D2(&mut self) -> f64 {
        self.xt2() * self.tT0D1().powi(2) + self.xt1() * self.tT0D2()
    }
    fn XT0D3(&mut self) -> f64 {
        self.xt3() * self.tT0D1().powi(3)
            + 3.0 * self.xt2() * self.tT0D1() * self.tT0D2()
            + self.xt1() * self.tT0D3()
    }
    fn XT0D4(&mut self) -> f64 {
        self.xt4() * self.tT0D1().powi(4)
            + 6.0 * self.xt3() * self.tT0D1().powi(2) * self.tT0D2()
            + 3.0 * self.xt2() * self.tT0D2().powi(2)
            + 4.0 * self.xt2() * self.tT0D1() * self.tT0D3()
            + self.xt1() * self.tT0D4()
    }
    fn XT1D0(&mut self) -> f64 {
        self.xt1() * self.tT1D0()
    }
    fn XT1D1(&mut self) -> f64 {
        self.xt2() * self.tT1D0() * self.tT0D1() + self.xt1() * self.tT1D1()
    }
    fn XT1D2(&mut self) -> f64 {
        self.xt3() * self.tT1D0() * self.tT0D1().powi(2)
            + 2.0 * self.xt2() * self.tT1D1() * self.tT0D1()
            + self.xt2() * self.tT1D0() * self.tT0D2()
            + self.xt1() * self.tT1D2()
    }
    fn XT1D3(&mut self) -> f64 {
        self.xt4() * self.tT1D0() * self.tT0D1().powi(3)
            + 3.0 * self.xt3() * self.tT1D1() * self.tT0D1().powi(2)
            + 3.0 * self.xt3() * self.tT1D0() * self.tT0D1() * self.tT0D2()
            + 3.0 * self.xt2() * self.tT1D2() * self.tT0D1()
            + 3.0 * self.xt2() * self.tT1D1() * self.tT0D2()
            + self.xt2() * self.tT1D0() * self.tT0D3()
            + self.xt1() * self.tT1D3()
    }
    fn XT2D0(&mut self) -> f64 {
        self.xt2() * self.tT1D0().powi(2) + self.xt1() * self.tT2D0()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn tT0D0(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * self.DeltaT0D0()
    }
    fn tT0D1(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * (self.DeltaT0D1() + self.DeltaT0D0())
    }
    fn tT0D2(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * (self.DeltaT0D2() + 2.0 * self.DeltaT0D1())
    }
    fn tT0D3(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * (self.DeltaT0D3() + 3.0 * self.DeltaT0D2())
    }
    fn tT0D4(&self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * (self.DeltaT0D4() + 4.0 * self.DeltaT0D3())
    }
    fn tT1D0(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * self.DeltaT1D0()
    }
    fn tT1D1(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * (self.DeltaT1D1() + self.DeltaT1D0())
    }
    fn tT1D2(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * (self.DeltaT1D2() + 2.0 * self.DeltaT1D1())
    }
    fn tT1D3(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * (self.DeltaT1D3() + 3.0 * self.DeltaT1D2())
    }
    fn tT2D0(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_sigma3 * self.DeltaT2D0()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn DeltaT0D0(&mut self) -> f64 {
        self.giiT0D0() * (self.epsilon_AB_temp.exp() - 1.0)
    }
    fn DeltaT0D1(&mut self) -> f64 {
        self.giiT0D1() * (self.epsilon_AB_temp.exp() - 1.0)
    }
    fn DeltaT0D2(&mut self) -> f64 {
        self.giiT0D2() * (self.epsilon_AB_temp.exp() - 1.0)
    }
    fn DeltaT0D3(&self) -> f64 {
        self.giiT0D3() * (self.epsilon_AB_temp.exp() - 1.0)
    }
    fn DeltaT0D4(&self) -> f64 {
        self.giiT0D4() * (self.epsilon_AB_temp.exp() - 1.0)
    }
    fn DeltaT1D0(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D0() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D0() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT1D1(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D1() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D1() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT1D2(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D2() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D2() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT1D3(&self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D3() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D3() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT2D0(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT2D0() * (epsilon_AB_T.exp() - 1.0)
            - 2.0 * self.giiT1D0() * epsilon_AB_T.exp() * epsilon_AB_T
            + self.giiT0D0() * epsilon_AB_T.exp() * epsilon_AB_T * (epsilon_AB_T + 2.0)
    }
}
const FRAC_RE30_NA: f64 = R / FRAC_NA_1E30;
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
