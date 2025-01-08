use thiserror::Error;
#[derive(Debug, Error)]
enum PcSaftPureErr {
    #[error("c_flash diverge")]
    NotConvForC,
    #[error("t_flash diverge")]
    NotConvForT,
    #[error("tp_flash diverge")]
    NotConvForTP,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
    #[error("property only in double phase")]
    OnlyInDoublePhase,
}
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_6, PI};
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
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)]
pub struct PcSaftPure {
    m: f64,
    sigma: f64,
    epsilon: f64,
    temp: f64,
    rho_num: f64,
    rhov_num: f64,
    rhol_num: f64,
    is_single_phase: bool,
    // changed from temp and rho_num
    eta: f64,
    eta1: f64,
    eta2: f64,
    m2e1s3: f64,
    m2e2s3: f64,
    // association term
    assoc_type: AssocType,
    epsilon_AB: f64,
    kappa_AB_plus: f64,
    // modified aly_lee_cv0
    mB: f64, // mB=B/R-1
    mC: f64, // mC=C/R
    mD: f64, // mD=D
    mE: f64, // mE=E/R
    mF: f64, // mF=F
    // A & B
    Ai0: f64,
    Ai1: f64,
    Ai2: f64,
    Ai3: f64,
    Ai4: f64,
    Ai5: f64,
    Ai6: f64,
    Bi0: f64,
    Bi1: f64,
    Bi2: f64,
    Bi3: f64,
    Bi4: f64,
    Bi5: f64,
    Bi6: f64,
    // cached variables
    giiT0D0: (f64, f64),
    giiT0D1: (f64, f64),
    giiT0D2: (f64, f64),
    cT0D0: (f64, f64),
    cT0D1: (f64, f64),
    cT0D2: (f64, f64),
    i2T0D0: (f64, f64),
    i2T0D1: (f64, f64),
    i2T0D2: (f64, f64),
}
impl PcSaftPure {
    pub fn new_fluid(m: f64, sigma: f64, epsilon: f64) -> Self {
        let m1 = (m - 1.0) / m; // (m-1)/m
        let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
        Self {
            m,
            sigma,
            epsilon,
            temp: 1.0,
            rho_num: 1E-10,
            rhov_num: 0.0,
            rhol_num: 0.0,
            is_single_phase: true,
            // changed from temp and rho_num
            eta: 0.0,
            eta1: 0.0,
            eta2: 0.0,
            m2e1s3: 0.0,
            m2e2s3: 0.0,
            // association term
            assoc_type: AssocType::Type0,
            epsilon_AB: 0.0,
            kappa_AB_plus: 0.0,
            // modified aly_lee_cv0
            mB: 1.5,
            mC: 0.0,
            mD: 1.0,
            mE: 0.0,
            mF: 1.0,
            // A & B
            Ai0: A00 + m1 * A10 + m12 * A20,
            Ai1: A01 + m1 * A11 + m12 * A21,
            Ai2: A02 + m1 * A12 + m12 * A22,
            Ai3: A03 + m1 * A13 + m12 * A23,
            Ai4: A04 + m1 * A14 + m12 * A24,
            Ai5: A05 + m1 * A15 + m12 * A25,
            Ai6: A06 + m1 * A16 + m12 * A26,
            Bi0: B00 + m1 * B10 + m12 * B20,
            Bi1: B01 + m1 * B11 + m12 * B21,
            Bi2: B02 + m1 * B12 + m12 * B22,
            Bi3: B03 + m1 * B13 + m12 * B23,
            Bi4: B04 + m1 * B14 + m12 * B24,
            Bi5: B05 + m1 * B15 + m12 * B25,
            Bi6: B06 + m1 * B16 + m12 * B26,
            // cached variables
            giiT0D0: (0.0, 0.0),
            giiT0D1: (0.0, 0.0),
            giiT0D2: (0.0, 0.0),
            cT0D0: (0.0, 0.0),
            cT0D1: (0.0, 0.0),
            cT0D2: (0.0, 0.0),
            i2T0D0: (0.0, 0.0),
            i2T0D1: (0.0, 0.0),
            i2T0D2: (0.0, 0.0),
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)]
impl PcSaftPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_1_assoc_type(&mut self, epsilon_AB: f64, kappa_AB: f64) {
        self.assoc_type = AssocType::Type1 { X: 1.0 };
        self.epsilon_AB = epsilon_AB;
        self.kappa_AB_plus = kappa_AB * self.sigma.powi(3);
    }
    pub fn set_2B_assoc_type(&mut self, epsilon_AB: f64, kappa_AB: f64) {
        self.assoc_type = AssocType::Type2B { X: 1.0 };
        self.epsilon_AB = epsilon_AB;
        self.kappa_AB_plus = kappa_AB * self.sigma.powi(3);
    }
    pub fn set_3B_assoc_type(&mut self, epsilon_AB: f64, kappa_AB: f64) {
        self.assoc_type = AssocType::Type3B { XA: 1.0 };
        self.epsilon_AB = epsilon_AB;
        self.kappa_AB_plus = kappa_AB * self.sigma.powi(3);
    }
    pub fn set_3Bm_assoc_type(&mut self, epsilon_AB: f64, kappa_AB: f64, b: f64) {
        self.assoc_type = AssocType::Type3Bm {
            XA: 1.0,
            b: b.abs().min(b.abs().recip()),
        };
        self.epsilon_AB = epsilon_AB;
        self.kappa_AB_plus = kappa_AB * self.sigma.powi(3);
    }
    pub fn set_aly_lee_cp0(&mut self, B: f64, C: f64, D: f64, E: f64, F: f64) {
        self.mB = B / R - 1.0; // mB = B/R -1
        self.mC = C / R; // mC = C/R
        self.mD = D; // mD = D
        self.mE = E / R; // mE = E/R
        self.mF = F; // mF = F
    }
    pub fn print_derivatives(&mut self) {
        self.check_derivatives(true);
    }
    pub fn td_unchecked(&mut self, temp: f64, rho_mol: f64) {
        self.set_temperature_and_number_density(temp, rho_mol * FRAC_NA_1E30);
        self.is_single_phase = true;
    }
    pub fn c_flash(&mut self) -> anyhow::Result<()> {
        // Iteration from temp_c = 1000 eta_c = 1E-10
        let mut temp_c = 1000.0;
        let mut dens_c = (1E-10)
            / (FRAC_PI_6
                * self.m
                * (self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / temp_c).exp())).powi(3));
        // Define variables
        let mut Dp_Drho_T = self.calc_Dp_Drho_T(temp_c, dens_c);
        let mut D2p_DTrho;
        let mut D2p_Drho2_T = self.calc_D2p_Drho2_T(temp_c, dens_c);
        let mut D3p_DTrho2;
        let mut D3p_Drho3_T;
        for _i in 1..100000 {
            D2p_DTrho = self.calc_D2p_DTrho(temp_c, dens_c);
            D3p_DTrho2 = self.calc_D3p_DTrho2(temp_c, dens_c);
            D3p_Drho3_T = self.calc_D3p_Drho3_T(temp_c, dens_c);
            temp_c -= (Dp_Drho_T * D3p_Drho3_T - D2p_Drho2_T * D2p_Drho2_T)
                / (D2p_DTrho * D3p_Drho3_T - D3p_DTrho2 * D2p_Drho2_T);
            dens_c -= (Dp_Drho_T * D3p_DTrho2 - D2p_Drho2_T * D2p_DTrho)
                / (D2p_Drho2_T * D3p_DTrho2 - D3p_Drho3_T * D2p_DTrho);
            Dp_Drho_T = self.calc_Dp_Drho_T(temp_c, dens_c);
            D2p_Drho2_T = self.calc_D2p_Drho2_T(temp_c, dens_c);
            if Dp_Drho_T.abs() < 1E3 && D2p_Drho2_T.abs() < 1E6 {
                self.set_temperature_and_number_density(temp_c, dens_c);
                self.is_single_phase = true;
                return Ok(());
            }
        }
        Err(anyhow!(PcSaftPureErr::NotConvForC))
    }
    pub fn t_flash(&mut self, temp: f64) -> anyhow::Result<()> {
        let d3 = (self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp())).powi(3);
        // Vapor phase: eta = 1E-10
        let rhov_num_guess = 1E-10 / (FRAC_PI_6 * self.m * d3);
        let mut rhov_num = rhov_num_guess;
        let mut Dp_Drhov_T = self.calc_Dp_Drho_T(temp, rhov_num);
        for _i in 1..100000 {
            if Dp_Drhov_T.abs() < 1.0 {
                break;
            } else {
                rhov_num -= Dp_Drhov_T / self.calc_D2p_Drho2_T(temp, rhov_num);
                Dp_Drhov_T = self.calc_Dp_Drho_T(temp, rhov_num);
            }
        }
        let pv_limit = self.calc_p(temp, rhov_num);
        if pv_limit.is_sign_negative() {
            return Err(anyhow!(PcSaftPureErr::NotConvForT));
        }
        // Liquid phase: eta = 0.5
        let rhol_num_guess = 0.5 / (FRAC_PI_6 * self.m * d3);
        let mut rhol_num = rhol_num_guess;
        let mut Dp_Drhol_T = self.calc_Dp_Drho_T(temp, rhol_num);
        for _i in 1..100000 {
            if Dp_Drhol_T.abs() < 1.0 {
                break;
            } else {
                rhol_num -= Dp_Drhol_T / self.calc_D2p_Drho2_T(temp, rhol_num);
                Dp_Drhol_T = self.calc_Dp_Drho_T(temp, rhol_num);
            }
        }
        if rhol_num < rhov_num {
            return Err(anyhow!(PcSaftPureErr::NotConvForT));
        }
        let pl_limit = self.calc_p(temp, rhol_num);
        if pl_limit > pv_limit {
            return Err(anyhow!(PcSaftPureErr::NotConvForT));
        }
        // Iteration for saturation state
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
            Err(anyhow!(PcSaftPureErr::NotConvForT))
        } else {
            self.is_single_phase = false;
            self.temp = temp;
            self.rhov_num = rhov_num;
            self.rhol_num = rhol_num;
            Ok(())
        }
    }
    pub fn tp_flash(&mut self, temp: f64, p: f64) -> anyhow::Result<()> {
        let d3 = (self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp())).powi(3);
        // Iteration from gas phase: eta = 1E-10
        let rhov_num = self.calc_density(temp, p, 1E-10 / (FRAC_PI_6 * self.m * d3));
        let lnphi_v = if rhov_num.is_nan() {
            f64::INFINITY
        } else {
            self.calc_lnphi(temp, rhov_num)
        };
        // Iteration from liquid phase: eta = 0.5
        let rhol_num = self.calc_density(temp, p, 0.5 / (FRAC_PI_6 * self.m * d3));
        let lnphi_l = if rhol_num.is_nan() {
            f64::INFINITY
        } else {
            self.calc_lnphi(temp, rhol_num)
        };
        // Select the correct output
        if lnphi_v.is_infinite() && lnphi_l.is_infinite() {
            Err(anyhow!(PcSaftPureErr::NotConvForTP))
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
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn rho(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.rho_num / FRAC_NA_1E30)
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn p(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_p(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn w(&mut self, M: f64) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok((self.calc_w2(self.temp, self.rho_num) / M).sqrt())
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn cv(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_cv(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn cp(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_cp(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn cp_res(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_cp_res(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn h_res(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_h_res(self.temp, self.rho_num))
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn p_s(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::OnlyInDoublePhase))
        } else {
            Ok(self.calc_p(self.temp, self.rhov_num) / 2.0
                + self.calc_p(self.temp, self.rhol_num) / 2.0)
        }
    }
    pub fn T_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::OnlyInDoublePhase))
        } else {
            Ok(self.temp)
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::OnlyInDoublePhase))
        } else {
            Ok(self.rhov_num / FRAC_NA_1E30)
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::OnlyInDoublePhase))
        } else {
            Ok(self.rhol_num / FRAC_NA_1E30)
        }
    }
    pub fn B(&mut self, temp: f64) -> f64 {
        let mut rho_num: f64 = 1E-9;
        let mut val_old: f64 = self.calc_rT0D1(temp, rho_num) / rho_num;
        let mut val_new: f64 = 0.0;
        loop {
            rho_num /= 10.0;
            if rho_num < 1E-30 {
                break;
            }
            val_new = self.calc_rT0D1(temp, rho_num) / rho_num;
            if (val_new / val_old - 1.0).abs() < 1E-9 {
                break;
            }
            val_old = val_new;
        }
        val_new * FRAC_NA_1E30
    }
    pub fn C(&mut self, temp: f64) -> f64 {
        let mut rho_num: f64 = 1E-9;
        let mut val_old: f64 = self.calc_rT0D2(temp, rho_num) / rho_num.powi(2);
        let mut val_new: f64 = 0.0;
        loop {
            rho_num /= 10.0;
            if rho_num < 1E-30 {
                break;
            }
            val_new = self.calc_rT0D2(temp, rho_num) / rho_num.powi(2);
            if (val_new / val_old - 1.0).abs() < 1E-9 {
                break;
            }
            val_old = val_new;
        }
        val_new * FRAC_NA_1E30.powi(2)
    }
    pub fn D(&mut self, temp: f64) -> f64 {
        let mut rho_num: f64 = 1E-9;
        let mut val_old: f64 = self.calc_rT0D3(temp, rho_num) / rho_num.powi(3);
        let mut val_new: f64 = 0.0;
        loop {
            rho_num /= 10.0;
            if rho_num < 1E-30 {
                break;
            }
            val_new = self.calc_rT0D3(temp, rho_num) / rho_num.powi(3);
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
        temp.into_iter()
            .zip(dens_num)
            .map(|(t, d)| self.calc_p(t, d))
            .collect()
    }
    pub fn vec_tp_flash(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
        temp.into_iter()
            .zip(pres)
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
        temp.into_iter()
            .zip(pres)
            .map(|(t, p)| {
                let d3 = (self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp())).powi(3);
                // Iteration from gas phase: eta = 1E-10
                self.calc_density(t, p, 1E-10 / (FRAC_PI_6 * self.m * d3))
            })
            .collect()
    }
    pub fn vec_tp_flash_l(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
        temp.into_iter()
            .zip(pres)
            .map(|(t, p)| {
                let d3 = (self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp())).powi(3);
                // Iteration from gas phase: eta = 0.5
                self.calc_density(t, p, 0.5 / (FRAC_PI_6 * self.m * d3))
            })
            .collect()
    }
    pub fn vec_cp(&mut self, temp: Vec<f64>, dens_num: Vec<f64>) -> Vec<f64> {
        temp.into_iter()
            .zip(dens_num)
            .map(|(t, d)| self.calc_cp(t, d))
            .collect()
    }
    pub fn vec_w(&mut self, temp: Vec<f64>, dens_num: Vec<f64>, M: f64) -> Vec<f64> {
        temp.into_iter()
            .zip(dens_num)
            .map(|(t, d)| (self.calc_w2(t, d) / M).sqrt()) // m/s
            .collect()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn calc_p(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp) * rho_num * (1.0 + self.calc_rT0D1(temp, rho_num))
    }
    fn calc_Dp_Drho_T(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (1.0 + 2.0 * self.calc_rT0D1(temp, rho_num) + self.calc_rT0D2(temp, rho_num))
    }
    fn calc_D2p_Drho2_T(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (2.0 * self.calc_rT0D1(temp, rho_num)
                + 4.0 * self.calc_rT0D2(temp, rho_num)
                + self.calc_rT0D3(temp, rho_num))
            / rho_num
    }
    fn calc_D3p_Drho3_T(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (6.0 * self.calc_rT0D2(temp, rho_num)
                + 6.0 * self.calc_rT0D3(temp, rho_num)
                + self.calc_rT0D4(temp, rho_num))
            / rho_num.powi(2)
    }
    fn calc_D2p_DTrho(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (1.0
                + 2.0 * self.calc_rT0D1(temp, rho_num)
                + self.calc_rT0D2(temp, rho_num)
                + 2.0 * self.calc_rT1D1(temp, rho_num)
                + self.calc_rT1D2(temp, rho_num))
    }
    fn calc_D3p_DTrho2(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (2.0 * self.calc_rT0D1(temp, rho_num)
                + 4.0 * self.calc_rT0D2(temp, rho_num)
                + self.calc_rT0D3(temp, rho_num)
                + 2.0 * self.calc_rT1D1(temp, rho_num)
                + 4.0 * self.calc_rT1D2(temp, rho_num)
                + self.calc_rT1D3(temp, rho_num))
            / rho_num
    }
    fn calc_w2(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * temp
            * ((1.0 + 2.0 * self.calc_rT0D1(temp, rho_num) + self.calc_rT0D2(temp, rho_num))
                + (1.0 + self.calc_rT0D1(temp, rho_num) + self.calc_rT1D1(temp, rho_num)).powi(2)
                    / (self.calc_ideal_cv(temp)
                        - 2.0 * self.calc_rT1D0(temp, rho_num)
                        - self.calc_rT2D0(temp, rho_num)))
    }
    fn calc_cv(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * (self.calc_ideal_cv(temp)
            - 2.0 * self.calc_rT1D0(temp, rho_num)
            - self.calc_rT2D0(temp, rho_num))
    }
    fn calc_cp(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * (self.calc_ideal_cv(temp)
            - 2.0 * self.calc_rT1D0(temp, rho_num)
            - self.calc_rT2D0(temp, rho_num)
            + (1.0 + self.calc_rT0D1(temp, rho_num) + self.calc_rT1D1(temp, rho_num)).powi(2)
                / (1.0 + 2.0 * self.calc_rT0D1(temp, rho_num) + self.calc_rT0D2(temp, rho_num)))
    }
    fn calc_cp_res(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * ((-1.0 - 2.0 * self.calc_rT1D0(temp, rho_num) - self.calc_rT2D0(temp, rho_num))
            + (1.0 + self.calc_rT0D1(temp, rho_num) + self.calc_rT1D1(temp, rho_num)).powi(2)
                / (1.0 + 2.0 * self.calc_rT0D1(temp, rho_num) + self.calc_rT0D2(temp, rho_num)))
    }
    fn calc_h_res(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * temp * (-self.calc_rT1D0(temp, rho_num) + self.calc_rT0D1(temp, rho_num))
    }
    fn calc_lnphi(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.calc_rT0D0(temp, rho_num) + self.calc_rT0D1(temp, rho_num)
            - (1.0 + self.calc_rT0D1(temp, rho_num)).ln()
    }
    fn calc_density(&mut self, temp: f64, p: f64, rho_num_guess: f64) -> f64 {
        let mut rho_num = rho_num_guess;
        let (mut p_diff, mut val_Dp_Drho_T, mut rho_num_diff);
        for _i in 1..10000 {
            p_diff = self.calc_p(temp, rho_num) - p;
            if p_diff.abs() < f64::EPSILON {
                return rho_num;
            }
            val_Dp_Drho_T = self.calc_Dp_Drho_T(temp, rho_num);
            if val_Dp_Drho_T.is_sign_negative() {
                return f64::NAN;
            }
            rho_num_diff = p_diff / val_Dp_Drho_T;
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
#[allow(non_snake_case)]
impl PcSaftPure {
    fn set_temperature_and_number_density(&mut self, temp: f64, rho_num: f64) {
        if temp != self.temp {
            self.temp = temp;
            self.m2e1s3 = self.m.powi(2) * (self.epsilon / temp) * self.sigma.powi(3);
            self.m2e2s3 = self.m.powi(2) * (self.epsilon / temp).powi(2) * self.sigma.powi(3);
        } else if rho_num != self.rho_num {
        } else {
            return;
        }
        self.rho_num = rho_num;
        let d = self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp());
        let d1 = -0.36 * self.sigma * (-3.0 * self.epsilon / temp).exp() * self.epsilon / temp;
        let d2 = d1 * (3.0 * self.epsilon / temp - 2.0);
        self.eta = FRAC_PI_6 * rho_num * self.m * d.powi(3);
        self.eta1 = FRAC_PI_2 * rho_num * self.m * d.powi(2) * d1;
        self.eta2 =
            PI * rho_num * self.m * d * d1.powi(2) + FRAC_PI_2 * rho_num * self.m * d.powi(2) * d2;
        if let AssocType::Type0 = self.assoc_type {
        } else {
            let t = self.tT0D0();
            match &mut self.assoc_type {
                AssocType::Type1 { X } | AssocType::Type2B { X } => {
                    *X = (-1.0 + (1.0 + 4.0 * t).sqrt()) / (2.0 * t);
                }
                AssocType::Type3B { XA } => {
                    *XA = (-(1.0 - t) + ((1.0 - t).powi(2) + 8.0 * t).sqrt()) / (4.0 * t);
                }
                AssocType::Type3Bm { XA, b } => {
                    let b = *b;
                    *XA = (-(1.0 - b * t) + ((1.0 - b * t).powi(2) + 4.0 * (1.0 + b) * t).sqrt())
                        / (2.0 * (1.0 + b) * t);
                }
                AssocType::Type0 => (),
            }
        }
    }
    fn calc_ideal_cv(&mut self, temp: f64) -> f64 {
        self.mB
            + self.mC * (self.mD / temp / (self.mD / temp).sinh()).powi(2)
            + self.mE * (self.mF / temp / (self.mF / temp).cosh()).powi(2)
    }
    fn calc_rT0D0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D0() + self.dispT0D0() + self.assocT0D0()
    }
    fn calc_rT0D1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D1() + self.dispT0D1() + self.assocT0D1()
    }
    fn calc_rT0D2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D2() + self.dispT0D2() + self.assocT0D2()
    }
    fn calc_rT0D3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D3() + self.dispT0D3() + self.assocT0D3()
    }
    fn calc_rT0D4(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT0D4() + self.dispT0D4() + self.assocT0D4()
    }
    fn calc_rT1D0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D0() + self.dispT1D0() + self.assocT1D0()
    }
    fn calc_rT1D1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D1() + self.dispT1D1() + self.assocT1D1()
    }
    fn calc_rT1D2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D2() + self.dispT1D2() + self.assocT1D2()
    }
    fn calc_rT1D3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hcT1D3() + self.dispT1D3() + self.assocT1D3()
    }
    fn calc_rT2D0(&mut self, temp: f64, rho_num: f64) -> f64 {
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
        let giiT0D2: f64 = self.giiT0D2();
        self.m * self.hsT0D4()
            + (1.0 - self.m)
                * (-6.0 / self.giiT0D0().powi(4) * self.giiT0D1().powi(4)
                    + 12.0 / self.giiT0D0().powi(3) * self.giiT0D1().powi(2) * giiT0D2
                    - 4.0 / self.giiT0D0().powi(2) * self.giiT0D1() * self.giiT0D3()
                    - 3.0 / self.giiT0D0().powi(2) * giiT0D2.powi(2)
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
        let giiT0D2: f64 = self.giiT0D2();
        let giiT1D0: f64 = self.giiT1D0();
        let giiT1D1: f64 = self.giiT1D1();
        self.m * self.hsT1D3()
            + (1.0 - self.m)
                * (-6.0 / self.giiT0D0().powi(4) * giiT1D0 * self.giiT0D1().powi(3)
                    + 6.0 / self.giiT0D0().powi(3) * giiT1D1 * self.giiT0D1().powi(2)
                    + 6.0 / self.giiT0D0().powi(3) * giiT1D0 * self.giiT0D1() * giiT0D2
                    - 3.0 / self.giiT0D0().powi(2) * giiT1D1 * giiT0D2
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
        self.eta1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
    }
    fn hsT1D1(&self) -> f64 {
        self.eta1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
    fn hsT1D2(&self) -> f64 {
        2.0 * self.eta1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
            + self.eta1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn hsT1D3(&self) -> f64 {
        3.0 * self.eta1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
            + self.eta1 * self.eta.powi(3) * (168.0 - 48.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn hsT2D0(&self) -> f64 {
        self.eta2 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta1.powi(2) * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn giiT0D0(&mut self) -> f64 {
        if self.eta != self.giiT0D0.0 {
            self.giiT0D0 = (self.eta, (1.0 - 0.5 * self.eta) / (1.0 - self.eta).powi(3))
        }
        self.giiT0D0.1
    }
    fn giiT0D1(&mut self) -> f64 {
        if self.eta != self.giiT0D1.0 {
            self.giiT0D1 = (
                self.eta,
                self.eta * (2.5 - self.eta) / (1.0 - self.eta).powi(4),
            )
        }
        self.giiT0D1.1
    }
    fn giiT0D2(&mut self) -> f64 {
        if self.eta != self.giiT0D2.0 {
            self.giiT0D2 = (
                self.eta,
                self.eta.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5),
            )
        }
        self.giiT0D2.1
    }
    fn giiT0D3(&self) -> f64 {
        self.eta.powi(3) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn giiT0D4(&self) -> f64 {
        self.eta.powi(4) * (240.0 - 60.0 * self.eta) / (1.0 - self.eta).powi(7)
    }
    fn giiT1D0(&self) -> f64 {
        self.eta1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
    }
    fn giiT1D1(&self) -> f64 {
        self.eta1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
            + self.eta1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn giiT1D2(&self) -> f64 {
        2.0 * self.eta1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
            + self.eta1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn giiT1D3(&self) -> f64 {
        3.0 * self.eta1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
            + self.eta1 * self.eta.powi(3) * (240.0 - 60.0 * self.eta) / (1.0 - self.eta).powi(7)
    }
    fn giiT2D0(&self) -> f64 {
        self.eta2 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
            + self.eta1.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
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
        let i2T0D3 = self.i2T0D3(self.eta);
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (4.0 * self.i1T0D3(self.eta) + self.i1T0D4(self.eta))
                + self.m
                    * self.m2e2s3
                    * (4.0 * c1T0D3 * self.i2T0D0()
                        + 12.0 * self.c1T0D2() * self.i2T0D1()
                        + 12.0 * self.c1T0D1() * self.i2T0D2()
                        + 4.0 * self.c1T0D0() * i2T0D3
                        + self.c1T0D4(self.eta) * self.i2T0D0()
                        + 4.0 * c1T0D3 * self.i2T0D1()
                        + 6.0 * self.c1T0D2() * self.i2T0D2()
                        + 4.0 * self.c1T0D1() * i2T0D3
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
        let i2T0D3 = self.i2T0D3(self.eta);
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
                        + c1T1D0 * i2T0D3
                        + self.c1T0D0() * self.i2T1D3(self.eta)
                        - 2.0 * self.c1T0D0() * i2T0D3))
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
        let cT0D2 = self.cT0D2();
        let cT1D0 = self.cT1D0(eta);
        let cT1D1 = self.cT1D1(eta);
        24.0 * cT1D0 * self.cT0D1().powi(3) / self.cT0D0().powi(5)
            - 18.0 * cT1D1 * self.cT0D1().powi(2) / self.cT0D0().powi(4)
            - 18.0 * cT1D0 * self.cT0D1() * cT0D2 / self.cT0D0().powi(4)
            + 6.0 * cT1D1 * cT0D2 / self.cT0D0().powi(3)
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
        if self.eta != self.cT0D0.0 {
            self.cT0D0 = (
                self.eta,
                1.0 + 2.0 * self.m * (4.0 * self.eta - self.eta.powi(2)) / (1.0 - self.eta).powi(4)
                    + (1.0 - self.m)
                        * (20.0 * self.eta - 27.0 * self.eta.powi(2) + 12.0 * self.eta.powi(3)
                            - 2.0 * self.eta.powi(4))
                        / ((1.0 - self.eta) * (2.0 - self.eta)).powi(2),
            )
        }
        self.cT0D0.1
    }
    fn cT0D1(&mut self) -> f64 {
        if self.eta != self.cT0D1.0 {
            self.cT0D1 = (
                self.eta,
                self.eta
                    * 2.0
                    * (2.0 * self.m * (-self.eta.powi(2) + 5.0 * self.eta + 2.0)
                        / (1.0 - self.eta).powi(5)
                        + (1.0 - self.m)
                            * (self.eta.powi(3) + 6.0 * self.eta.powi(2) - 24.0 * self.eta + 20.0)
                            / ((1.0 - self.eta) * (2.0 - self.eta)).powi(3)),
            )
        }
        self.cT0D1.1
    }
    fn cT0D2(&mut self) -> f64 {
        if self.eta != self.cT0D2.0 {
            self.cT0D2 = (
                self.eta,
                self.eta.powi(2)
                    * 6.0
                    * (2.0 * self.m * (-self.eta.powi(2) + 6.0 * self.eta + 5.0)
                        / (1.0 - self.eta).powi(6)
                        + (1.0 - self.m)
                            * (-self.eta.powi(4) - 8.0 * self.eta.powi(3)
                                + 48.0 * self.eta.powi(2)
                                - 80.0 * self.eta
                                + 44.0)
                            / ((1.0 - self.eta) * (2.0 - self.eta)).powi(4)),
            )
        }
        self.cT0D2.1
    }
    fn cT0D3(&self, eta: f64) -> f64 {
        eta.powi(3)
            * 24.0
            * (2.0 * self.m * (-eta.powi(2) + 7.0 * eta + 9.0) / (1.0 - eta).powi(7)
                + (1.0 - self.m)
                    * (eta.powi(5) + 10.0 * eta.powi(4) - 80.0 * eta.powi(3) + 200.0 * eta.powi(2)
                        - 220.0 * eta
                        + 92.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(5))
    }
    fn cT0D4(&self, eta: f64) -> f64 {
        eta.powi(4)
            * 120.0
            * (2.0 * self.m * (-eta.powi(2) + 8.0 * eta + 14.0) / (1.0 - eta).powi(8)
                + (1.0 - self.m)
                    * (-eta.powi(6) - 12.0 * eta.powi(5) + 120.0 * eta.powi(4)
                        - 400.0 * eta.powi(3)
                        + 660.0 * eta.powi(2)
                        - 552.0 * eta
                        + 188.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(6))
    }
    fn cT1D0(&self, eta: f64) -> f64 {
        self.eta1
            * 2.0
            * (2.0 * self.m * (-eta.powi(2) + 5.0 * eta + 2.0) / (1.0 - eta).powi(5)
                + (1.0 - self.m) * (eta.powi(3) + 6.0 * eta.powi(2) - 24.0 * eta + 20.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(3))
    }
    fn cT1D1(&self, eta: f64) -> f64 {
        self.eta1
            * 2.0
            * (2.0 * self.m * (-eta.powi(2) + 5.0 * eta + 2.0) / (1.0 - eta).powi(5)
                + (1.0 - self.m) * (eta.powi(3) + 6.0 * eta.powi(2) - 24.0 * eta + 20.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(3))
            + self.eta1
                * self.eta
                * 6.0
                * (2.0 * self.m * (-eta.powi(2) + 6.0 * eta + 5.0) / (1.0 - eta).powi(6)
                    + (1.0 - self.m)
                        * (-eta.powi(4) - 8.0 * eta.powi(3) + 48.0 * eta.powi(2) - 80.0 * eta
                            + 44.0)
                        / ((1.0 - eta) * (2.0 - eta)).powi(4))
    }
    fn cT1D2(&self, eta: f64) -> f64 {
        self.eta1
            * self.eta
            * 12.0
            * (2.0 * self.m * (-eta.powi(2) + 6.0 * eta + 5.0) / (1.0 - eta).powi(6)
                + (1.0 - self.m)
                    * (-eta.powi(4) - 8.0 * eta.powi(3) + 48.0 * eta.powi(2) - 80.0 * eta + 44.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(4))
            + self.eta1
                * self.eta.powi(2)
                * 24.0
                * (2.0 * self.m * (-eta.powi(2) + 7.0 * eta + 9.0) / (1.0 - eta).powi(7)
                    + (1.0 - self.m)
                        * (eta.powi(5) + 10.0 * eta.powi(4) - 80.0 * eta.powi(3)
                            + 200.0 * eta.powi(2)
                            - 220.0 * eta
                            + 92.0)
                        / ((1.0 - eta) * (2.0 - eta)).powi(5))
    }
    fn cT1D3(&self, eta: f64) -> f64 {
        self.eta1
            * self.eta.powi(2)
            * 72.0
            * (2.0 * self.m * (-eta.powi(2) + 7.0 * eta + 9.0) / (1.0 - eta).powi(7)
                + (1.0 - self.m)
                    * (eta.powi(5) + 10.0 * eta.powi(4) - 80.0 * eta.powi(3) + 200.0 * eta.powi(2)
                        - 220.0 * eta
                        + 92.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(5))
            + self.eta1
                * self.eta.powi(3)
                * 120.0
                * (2.0 * self.m * (-eta.powi(2) + 8.0 * eta + 14.0) / (1.0 - eta).powi(8)
                    + (1.0 - self.m)
                        * (-eta.powi(6) - 12.0 * eta.powi(5) + 120.0 * eta.powi(4)
                            - 400.0 * eta.powi(3)
                            + 660.0 * eta.powi(2)
                            - 552.0 * eta
                            + 188.0)
                        / ((1.0 - eta) * (2.0 - eta)).powi(6))
    }
    fn cT2D0(&self, eta: f64) -> f64 {
        self.eta2
            * 2.0
            * (2.0 * self.m * (-eta.powi(2) + 5.0 * eta + 2.0) / (1.0 - eta).powi(5)
                + (1.0 - self.m) * (eta.powi(3) + 6.0 * eta.powi(2) - 24.0 * eta + 20.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(3))
            + self.eta1.powi(2)
                * 6.0
                * (2.0 * self.m * (-eta.powi(2) + 6.0 * eta + 5.0) / (1.0 - eta).powi(6)
                    + (1.0 - self.m)
                        * (-eta.powi(4) - 8.0 * eta.powi(3) + 48.0 * eta.powi(2) - 80.0 * eta
                            + 44.0)
                        / ((1.0 - eta) * (2.0 - eta)).powi(4))
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn i1T0D0(&self, eta: f64) -> f64 {
        self.Ai0
            + self.Ai1 * eta
            + self.Ai2 * eta.powi(2)
            + self.Ai3 * eta.powi(3)
            + self.Ai4 * eta.powi(4)
            + self.Ai5 * eta.powi(5)
            + self.Ai6 * eta.powi(6)
    }
    fn i1T0D1(&self, eta: f64) -> f64 {
        self.Ai1 * eta
            + self.Ai2 * 2.0 * eta.powi(2)
            + self.Ai3 * 3.0 * eta.powi(3)
            + self.Ai4 * 4.0 * eta.powi(4)
            + self.Ai5 * 5.0 * eta.powi(5)
            + self.Ai6 * 6.0 * eta.powi(6)
    }
    fn i1T0D2(&self, eta: f64) -> f64 {
        self.Ai2 * 2.0 * eta.powi(2)
            + self.Ai3 * 6.0 * eta.powi(3)
            + self.Ai4 * 12.0 * eta.powi(4)
            + self.Ai5 * 20.0 * eta.powi(5)
            + self.Ai6 * 30.0 * eta.powi(6)
    }
    fn i1T0D3(&self, eta: f64) -> f64 {
        self.Ai3 * 6.0 * eta.powi(3)
            + self.Ai4 * 24.0 * eta.powi(4)
            + self.Ai5 * 60.0 * eta.powi(5)
            + self.Ai6 * 120.0 * eta.powi(6)
    }
    fn i1T0D4(&self, eta: f64) -> f64 {
        self.Ai4 * 24.0 * eta.powi(4)
            + self.Ai5 * 120.0 * eta.powi(5)
            + self.Ai6 * 360.0 * eta.powi(6)
    }
    fn i1T1D0(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Ai1
                + self.Ai2 * 2.0 * eta
                + self.Ai3 * 3.0 * eta.powi(2)
                + self.Ai4 * 4.0 * eta.powi(3)
                + self.Ai5 * 5.0 * eta.powi(4)
                + self.Ai6 * 6.0 * eta.powi(5))
    }
    fn i1T1D1(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Ai1
                + self.Ai2 * 4.0 * eta
                + self.Ai3 * 9.0 * eta.powi(2)
                + self.Ai4 * 16.0 * eta.powi(3)
                + self.Ai5 * 25.0 * eta.powi(4)
                + self.Ai6 * 36.0 * eta.powi(5))
    }
    fn i1T1D2(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Ai2 * 4.0 * eta
                + self.Ai3 * 18.0 * eta.powi(2)
                + self.Ai4 * 48.0 * eta.powi(3)
                + self.Ai5 * 100.0 * eta.powi(4)
                + self.Ai6 * 180.0 * eta.powi(5))
    }
    fn i1T1D3(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Ai3 * 18.0 * eta.powi(2)
                + self.Ai4 * 96.0 * eta.powi(3)
                + self.Ai5 * 300.0 * eta.powi(4)
                + self.Ai6 * 720.0 * eta.powi(5))
    }
    fn i1T2D0(&self, eta: f64) -> f64 {
        self.eta1.powi(2)
            * (self.Ai2 * 2.0
                + self.Ai3 * 6.0 * eta
                + self.Ai4 * 12.0 * eta.powi(2)
                + self.Ai5 * 20.0 * eta.powi(3)
                + self.Ai6 * 30.0 * eta.powi(4))
            + self.eta2
                * (self.Ai1
                    + self.Ai2 * 2.0 * eta
                    + self.Ai3 * 3.0 * eta.powi(2)
                    + self.Ai4 * 4.0 * eta.powi(3)
                    + self.Ai5 * 5.0 * eta.powi(4)
                    + self.Ai6 * 6.0 * eta.powi(5))
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn i2T0D0(&mut self) -> f64 {
        if self.eta != self.i2T0D0.0 {
            self.i2T0D0 = (
                self.eta,
                self.Bi0
                    + self.Bi1 * self.eta
                    + self.Bi2 * self.eta.powi(2)
                    + self.Bi3 * self.eta.powi(3)
                    + self.Bi4 * self.eta.powi(4)
                    + self.Bi5 * self.eta.powi(5)
                    + self.Bi6 * self.eta.powi(6),
            )
        }
        self.i2T0D0.1
    }
    fn i2T0D1(&mut self) -> f64 {
        if self.eta != self.i2T0D1.0 {
            self.i2T0D1 = (
                self.eta,
                self.Bi1 * self.eta
                    + self.Bi2 * 2.0 * self.eta.powi(2)
                    + self.Bi3 * 3.0 * self.eta.powi(3)
                    + self.Bi4 * 4.0 * self.eta.powi(4)
                    + self.Bi5 * 5.0 * self.eta.powi(5)
                    + self.Bi6 * 6.0 * self.eta.powi(6),
            )
        }
        self.i2T0D1.1
    }
    fn i2T0D2(&mut self) -> f64 {
        if self.eta != self.i2T0D2.0 {
            self.i2T0D2 = (
                self.eta,
                self.Bi2 * 2.0 * self.eta.powi(2)
                    + self.Bi3 * 6.0 * self.eta.powi(3)
                    + self.Bi4 * 12.0 * self.eta.powi(4)
                    + self.Bi5 * 20.0 * self.eta.powi(5)
                    + self.Bi6 * 30.0 * self.eta.powi(6),
            )
        }
        self.i2T0D2.1
    }
    fn i2T0D3(&self, eta: f64) -> f64 {
        self.Bi3 * 6.0 * eta.powi(3)
            + self.Bi4 * 24.0 * eta.powi(4)
            + self.Bi5 * 60.0 * eta.powi(5)
            + self.Bi6 * 120.0 * eta.powi(6)
    }
    fn i2T0D4(&self, eta: f64) -> f64 {
        self.Bi4 * 24.0 * eta.powi(4)
            + self.Bi5 * 120.0 * eta.powi(5)
            + self.Bi6 * 360.0 * eta.powi(6)
    }
    fn i2T1D0(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Bi1
                + self.Bi2 * 2.0 * eta
                + self.Bi3 * 3.0 * eta.powi(2)
                + self.Bi4 * 4.0 * eta.powi(3)
                + self.Bi5 * 5.0 * eta.powi(4)
                + self.Bi6 * 6.0 * eta.powi(5))
    }
    fn i2T1D1(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Bi1
                + self.Bi2 * 4.0 * eta
                + self.Bi3 * 9.0 * eta.powi(2)
                + self.Bi4 * 16.0 * eta.powi(3)
                + self.Bi5 * 25.0 * eta.powi(4)
                + self.Bi6 * 36.0 * eta.powi(5))
    }
    fn i2T1D2(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Bi2 * 4.0 * eta
                + self.Bi3 * 18.0 * eta.powi(2)
                + self.Bi4 * 48.0 * eta.powi(3)
                + self.Bi5 * 100.0 * eta.powi(4)
                + self.Bi6 * 180.0 * eta.powi(5))
    }
    fn i2T1D3(&self, eta: f64) -> f64 {
        self.eta1
            * (self.Bi3 * 18.0 * eta.powi(2)
                + self.Bi4 * 96.0 * eta.powi(3)
                + self.Bi5 * 300.0 * eta.powi(4)
                + self.Bi6 * 720.0 * eta.powi(5))
    }
    fn i2T2D0(&self, eta: f64) -> f64 {
        self.eta1.powi(2)
            * (self.Bi2 * 2.0
                + self.Bi3 * 6.0 * eta
                + self.Bi4 * 12.0 * eta.powi(2)
                + self.Bi5 * 20.0 * eta.powi(3)
                + self.Bi6 * 30.0 * eta.powi(4))
            + self.eta2
                * (self.Bi1
                    + self.Bi2 * 2.0 * eta
                    + self.Bi3 * 3.0 * eta.powi(2)
                    + self.Bi4 * 4.0 * eta.powi(3)
                    + self.Bi5 * 5.0 * eta.powi(4)
                    + self.Bi6 * 6.0 * eta.powi(5))
    }
}
#[allow(non_snake_case)]
enum AssocType {
    Type0,
    Type1 { X: f64 },
    Type2B { X: f64 },
    Type3B { XA: f64 },
    Type3Bm { XA: f64, b: f64 },
}
impl AssocType {
    fn t1(&self) -> f64 {
        match self {
            AssocType::Type1 { X } | AssocType::Type2B { X } => X.powi(3) / (X - 2.0),
            AssocType::Type3B { XA } => {
                (XA * (2.0 * XA - 1.0)).powi(2) / (2.0 * XA.powi(2) - 4.0 * XA + 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                (XA * (XA + b * XA - b)).powi(2)
                    / ((1.0 + b) * XA.powi(2) - 2.0 * (1.0 + b) * XA + b)
            }
            _ => 0.0,
        }
    }
    fn t2(&self) -> f64 {
        2.0 * match self {
            AssocType::Type1 { X } | AssocType::Type2B { X } => {
                X.powi(5) / (X - 2.0).powi(3) * (X - 3.0)
            }
            AssocType::Type3B { XA } => {
                (XA * (2.0 * XA - 1.0)).powi(3) / (2.0 * XA.powi(2) - 4.0 * XA + 1.0).powi(3)
                    * (4.0 * XA.powi(3) - 12.0 * XA.powi(2) + 6.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                (XA * (XA + b * XA - b)).powi(3)
                    / ((1.0 + b) * XA.powi(2) - 2.0 * (1.0 + b) * XA + b).powi(3)
                    * ((1.0 + b).powi(2) * XA.powi(3) - 3.0 * (1.0 + b).powi(2) * XA.powi(2)
                        + (3.0 * b * (1.0 + b) * XA - b.powi(2)))
            }
            _ => 0.0,
        }
    }
    fn t3(&self) -> f64 {
        6.0 * match self {
            AssocType::Type1 { X } | AssocType::Type2B { X } => {
                X.powi(7) / (X - 2.0).powi(5) * (X.powi(2) - 6.0 * X + 10.0)
            }
            AssocType::Type3B { XA } => {
                (XA * (2.0 * XA - 1.0)).powi(4) / (2.0 * XA.powi(2) - 4.0 * XA + 1.0).powi(5)
                    * (16.0 * XA.powi(6) - 96.0 * XA.powi(5)
                        + (200.0 * XA.powi(4) - 160.0 * XA.powi(3))
                        + (62.0 * XA.powi(2) - 12.0 * XA + 1.0))
            }
            AssocType::Type3Bm { XA, b } => {
                (XA * (XA + b * XA - b)).powi(4)
                    / ((1.0 + b) * XA.powi(2) - 2.0 * (1.0 + b) * XA + b).powi(5)
                    * ((1.0 + b).powi(4) * XA.powi(6) - 6.0 * (1.0 + b).powi(4) * XA.powi(5)
                        + 5.0 * (3.0 * b + 2.0) * (1.0 + b).powi(3) * XA.powi(4)
                        - 20.0 * b * (1.0 + b).powi(3) * XA.powi(3)
                        + b.powi(2) * (15.0 * b + 16.0) * (1.0 + b) * XA.powi(2)
                        - (6.0 * b.powi(3) * (1.0 + b) * XA - b.powi(4)))
            }
            _ => 0.0,
        }
    }
    fn t4(&self) -> f64 {
        24.0 * match self {
            AssocType::Type1 { X } | AssocType::Type2B { X } => {
                X.powi(9) / (X - 2.0).powi(7) * (X.powi(3) - 9.0 * X.powi(2) + 29.0 * X - 35.0)
            }
            AssocType::Type3B { XA } => {
                (XA * (2.0 * XA - 1.0)).powi(5) / (2.0 * XA.powi(2) - 4.0 * XA + 1.0).powi(7)
                    * (64.0 * XA.powi(9) - 576.0 * XA.powi(8)
                        + (2080.0 * XA.powi(7) - 3808.0 * XA.powi(6))
                        + (3696.0 * XA.powi(5) - 2084.0 * XA.powi(4))
                        + (716.0 * XA.powi(3) - 150.0 * XA.powi(2))
                        + (18.0 * XA - 1.0))
            }
            AssocType::Type3Bm { XA, b } => {
                (XA * (XA + b * XA - b)).powi(5)
                    / ((1.0 + b) * XA.powi(2) - 2.0 * (1.0 + b) * XA + b).powi(7)
                    * ((1.0 + b).powi(6) * XA.powi(9) - 9.0 * (1.0 + b).powi(6) * XA.powi(8)
                        + (36.0 * b + 29.0) * (1.0 + b).powi(5) * XA.powi(7)
                        - 7.0 * (12.0 * b + 5.0) * (1.0 + b).powi(5) * XA.powi(6)
                        + 21.0 * b * (6.0 * b + 5.0) * (1.0 + b).powi(4) * XA.powi(5)
                        - b.powi(2)
                            * (126.0 * b.powi(2) + 260.0 * b + 135.0)
                            * (1.0 + b).powi(2)
                            * XA.powi(4)
                        + b.powi(3) * (84.0 * b + 95.0) * (1.0 + b).powi(2) * XA.powi(3)
                        - 3.0 * b.powi(4) * (12.0 * b + 13.0) * (1.0 + b) * XA.powi(2)
                        + (9.0 * b.powi(5) * (1.0 + b) * XA - b.powi(6)))
            }
            _ => 0.0,
        }
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn assocT0D0(&self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => X.ln() - X / 2.0 + 0.5,
            AssocType::Type2B { X } => 2.0 * X.ln() - X + 1.0,
            AssocType::Type3B { XA } => 2.0 * XA.ln() + (2.0 * XA - 1.0).ln() - 2.0 * XA + 2.0,
            AssocType::Type3Bm { XA, b } => {
                XA.ln() + (b * XA + 1.0 - b).ln() + ((1.0 + b) * XA - b).ln() - (1.0 + b) * XA
                    + (1.0 + b)
            }
        }
    }
    fn assocT0D1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT0D1(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT0D1(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT0D1(1.0, XA) + self.siteT0D1(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT0D1(1.0, XA)
                    + self.siteT0D1(b, b * XA + 1.0 - b)
                    + self.siteT0D1(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT0D2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT0D2(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT0D2(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT0D2(1.0, XA) + self.siteT0D2(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT0D2(1.0, XA)
                    + self.siteT0D2(b, b * XA + 1.0 - b)
                    + self.siteT0D2(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT0D3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT0D3(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT0D3(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT0D3(1.0, XA) + self.siteT0D3(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT0D3(1.0, XA)
                    + self.siteT0D3(b, b * XA + 1.0 - b)
                    + self.siteT0D3(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT0D4(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT0D4(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT0D4(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT0D4(1.0, XA) + self.siteT0D4(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT0D4(1.0, XA)
                    + self.siteT0D4(b, b * XA + 1.0 - b)
                    + self.siteT0D4(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT1D0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT1D0(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT1D0(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT1D0(1.0, XA) + self.siteT1D0(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT1D0(1.0, XA)
                    + self.siteT1D0(b, b * XA + 1.0 - b)
                    + self.siteT1D0(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT1D1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT1D1(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT1D1(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT1D1(1.0, XA) + self.siteT1D1(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT1D1(1.0, XA)
                    + self.siteT1D1(b, b * XA + 1.0 - b)
                    + self.siteT1D1(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT1D2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT1D2(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT1D2(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT1D2(1.0, XA) + self.siteT1D2(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT1D2(1.0, XA)
                    + self.siteT1D2(b, b * XA + 1.0 - b)
                    + self.siteT1D2(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT1D3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT1D3(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT1D3(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT1D3(1.0, XA) + self.siteT1D3(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT1D3(1.0, XA)
                    + self.siteT1D3(b, b * XA + 1.0 - b)
                    + self.siteT1D3(1.0 + b, (1.0 + b) * XA - b)
            }
        }
    }
    fn assocT2D0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => self.siteT2D0(1.0, X),
            AssocType::Type2B { X } => 2.0 * self.siteT2D0(1.0, X),
            AssocType::Type3B { XA } => {
                2.0 * self.siteT2D0(1.0, XA) + self.siteT2D0(2.0, 2.0 * XA - 1.0)
            }
            AssocType::Type3Bm { XA, b } => {
                self.siteT2D0(1.0, XA)
                    + self.siteT2D0(b, b * XA + 1.0 - b)
                    + self.siteT2D0(1.0 + b, (1.0 + b) * XA - b)
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
#[allow(non_snake_case)]
impl PcSaftPure {
    fn XT0D1(&mut self) -> f64 {
        self.assoc_type.t1() * self.tT0D1()
    }
    fn XT0D2(&mut self) -> f64 {
        self.assoc_type.t2() * self.tT0D1().powi(2) + self.assoc_type.t1() * self.tT0D2()
    }
    fn XT0D3(&mut self) -> f64 {
        self.assoc_type.t3() * self.tT0D1().powi(3)
            + 3.0 * self.assoc_type.t2() * self.tT0D1() * self.tT0D2()
            + self.assoc_type.t1() * self.tT0D3()
    }
    fn XT0D4(&mut self) -> f64 {
        self.assoc_type.t4() * self.tT0D1().powi(4)
            + 6.0 * self.assoc_type.t3() * self.tT0D1().powi(2) * self.tT0D2()
            + 3.0 * self.assoc_type.t2() * self.tT0D2().powi(2)
            + 4.0 * self.assoc_type.t2() * self.tT0D1() * self.tT0D3()
            + self.assoc_type.t1() * self.tT0D4()
    }
    fn XT1D0(&mut self) -> f64 {
        self.assoc_type.t1() * self.tT1D0()
    }
    fn XT1D1(&mut self) -> f64 {
        self.assoc_type.t2() * self.tT1D0() * self.tT0D1() + self.assoc_type.t1() * self.tT1D1()
    }
    fn XT1D2(&mut self) -> f64 {
        self.assoc_type.t3() * self.tT1D0() * self.tT0D1().powi(2)
            + 2.0 * self.assoc_type.t2() * self.tT1D1() * self.tT0D1()
            + self.assoc_type.t2() * self.tT1D0() * self.tT0D2()
            + self.assoc_type.t1() * self.tT1D2()
    }
    fn XT1D3(&mut self) -> f64 {
        self.assoc_type.t4() * self.tT1D0() * self.tT0D1().powi(3)
            + 3.0 * self.assoc_type.t3() * self.tT1D1() * self.tT0D1().powi(2)
            + 3.0 * self.assoc_type.t3() * self.tT1D0() * self.tT0D1() * self.tT0D2()
            + 3.0 * self.assoc_type.t2() * self.tT1D2() * self.tT0D1()
            + 3.0 * self.assoc_type.t2() * self.tT1D1() * self.tT0D2()
            + self.assoc_type.t2() * self.tT1D0() * self.tT0D3()
            + self.assoc_type.t1() * self.tT1D3()
    }
    fn XT2D0(&mut self) -> f64 {
        self.assoc_type.t2() * self.tT1D0().powi(2) + self.assoc_type.t1() * self.tT2D0()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn tT0D0(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * self.DeltaT0D0()
    }
    fn tT0D1(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * (self.DeltaT0D1() + self.DeltaT0D0())
    }
    fn tT0D2(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * (self.DeltaT0D2() + 2.0 * self.DeltaT0D1())
    }
    fn tT0D3(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * (self.DeltaT0D3() + 3.0 * self.DeltaT0D2())
    }
    fn tT0D4(&self) -> f64 {
        self.rho_num * self.kappa_AB_plus * (self.DeltaT0D4() + 4.0 * self.DeltaT0D3())
    }
    fn tT1D0(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * self.DeltaT1D0()
    }
    fn tT1D1(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * (self.DeltaT1D1() + self.DeltaT1D0())
    }
    fn tT1D2(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * (self.DeltaT1D2() + 2.0 * self.DeltaT1D1())
    }
    fn tT1D3(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * (self.DeltaT1D3() + 3.0 * self.DeltaT1D2())
    }
    fn tT2D0(&mut self) -> f64 {
        self.rho_num * self.kappa_AB_plus * self.DeltaT2D0()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn DeltaT0D0(&mut self) -> f64 {
        self.giiT0D0() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D1(&mut self) -> f64 {
        self.giiT0D1() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D2(&mut self) -> f64 {
        self.giiT0D2() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D3(&self) -> f64 {
        self.giiT0D3() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D4(&self) -> f64 {
        self.giiT0D4() * ((self.epsilon_AB / self.temp).exp() - 1.0)
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
const FRAC_NA_1E30: f64 = 6.02214076E-7; // const NA: f64 = 6.02214076E23;
const R: f64 = 8.314462618;
const A00: f64 = 0.9105631445;
const A01: f64 = 0.6361281449;
const A02: f64 = 2.6861347891;
const A03: f64 = -26.547362491;
const A04: f64 = 97.759208784;
const A05: f64 = -159.59154087;
const A06: f64 = 91.297774084;
const A10: f64 = -0.3084016918;
const A11: f64 = 0.1860531159;
const A12: f64 = -2.5030047259;
const A13: f64 = 21.419793629;
const A14: f64 = -65.255885330;
const A15: f64 = 83.318680481;
const A16: f64 = -33.746922930;
const A20: f64 = -0.0906148351;
const A21: f64 = 0.4527842806;
const A22: f64 = 0.5962700728;
const A23: f64 = -1.7241829131;
const A24: f64 = -4.1302112531;
const A25: f64 = 13.776631870;
const A26: f64 = -8.6728470368;
const B00: f64 = 0.7240946941;
const B01: f64 = 2.2382791861;
const B02: f64 = -4.0025849485;
const B03: f64 = -21.003576815;
const B04: f64 = 26.855641363;
const B05: f64 = 206.55133841;
const B06: f64 = -355.60235612;
const B10: f64 = -0.5755498075;
const B11: f64 = 0.6995095521;
const B12: f64 = 3.8925673390;
const B13: f64 = -17.215471648;
const B14: f64 = 192.67226447;
const B15: f64 = -161.82646165;
const B16: f64 = -165.20769346;
const B20: f64 = 0.0976883116;
const B21: f64 = -0.2557574982;
const B22: f64 = -9.1558561530;
const B23: f64 = 20.642075974;
const B24: f64 = -38.804430052;
const B25: f64 = 93.626774077;
const B26: f64 = -29.666905585;
impl PcSaftPure {
    pub fn check_derivatives(&mut self, print_val: bool) {
        let temp = self.temp;
        let rho = self.rho_num;
        if print_val {
            println!(
                "[rT0D0 == rT0D0] calc_rT0D0() ={}",
                self.calc_rT0D0(temp, rho)
            );
        }
        let compare_val = |val_calc: f64, val_diff: f64| {
            assert_eq!(
                &val_calc.abs().to_string()[0..11],
                &val_diff.abs().to_string()[0..11]
            )
        };
        // derivative for density
        let val_calc = self.calc_rT0D1(temp, rho) / rho;
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT0D0(temp, rhox), rho);
        if print_val {
            println!("[rT0D1 == rT0D1] calc_rT0D1() ={}", val_calc);
            println!("[rT0D0 -> rT0D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density
        let val_calc = self.calc_rT0D2(temp, rho) / rho.powi(2);
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT0D1(temp, rhox) / rhox, rho);
        if print_val {
            println!("[rT0D2 == rT0D2] calc_rT0D2() ={}", val_calc);
            println!("[rT0D1 -> rT0D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density
        let val_calc = self.calc_rT0D3(temp, rho) / rho.powi(3);
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT0D2(temp, rhox) / rhox.powi(2), rho);
        if print_val {
            println!("[rT0D3 == rT0D3] calc_rT0D3() ={}", val_calc);
            println!("[rT0D2 -> rT0D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density+density
        let val_calc = self.calc_rT0D4(temp, rho) / rho.powi(4);
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT0D3(temp, rhox) / rhox.powi(3), rho);
        if print_val {
            println!("[rT0D4 == rT0D4] calc_rT0D4() ={}", val_calc);
            println!("[rT0D3 -> rT0D4] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature
        let val_calc = self.calc_rT1D0(temp, rho) / temp;
        let val_diff = romberg_diff(|tempx: f64| self.calc_rT0D0(tempx, rho), temp);
        if print_val {
            println!("[rT1D0 == rT1D0] calc_rT1D0() ={}", val_calc);
            println!("[rT0D0 -> rT1D0] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density
        let val_calc = self.calc_rT1D1(temp, rho) / temp / rho;
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT1D0(temp, rhox) / temp, rho);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_rT1D1() ={}", val_calc);
            println!("[rT1D0 -> rT1D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tempx: f64| self.calc_rT0D1(tempx, rho) / rho, temp);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_rT1D1() ={}", val_calc);
            println!("[rT0D1 -> rT1D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density
        let val_calc = self.calc_rT1D2(temp, rho) / temp / rho.powi(2);
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT1D1(temp, rhox) / temp / rhox, rho);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_rT1D2() ={}", val_calc);
            println!("[rT1D1 -> rT1D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tempx: f64| self.calc_rT0D2(tempx, rho) / rho.powi(2), temp);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_rT1D2() ={}", val_calc);
            println!("[rT0D2 -> rT1D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density+density
        let val_calc = self.calc_rT1D3(temp, rho) / temp / rho.powi(3);
        let val_diff = romberg_diff(
            |rhox: f64| self.calc_rT1D2(temp, rhox) / temp / rhox.powi(2),
            rho,
        );
        if print_val {
            println!("[rT1D3 == rT1D3] calc_rT1D3() ={}", val_calc);
            println!("[rT1D2 -> rT1D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tempx: f64| self.calc_rT0D3(tempx, rho) / rho.powi(3), temp);
        if print_val {
            println!("[rT1D3 == rT1D3] calc_rT1D3() ={}", val_calc);
            println!("[rT0D3 -> rT1D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+temperature
        let val_calc = self.calc_rT2D0(temp, rho) / temp.powi(2);
        let val_diff = romberg_diff(|tempx: f64| self.calc_rT1D0(tempx, rho) / tempx, temp);
        if print_val {
            println!("[rT2D0 == rT2D0] calc_rT2D0() ={}", val_calc);
            println!("[rT1D0 -> rT2D0] romberg_diff ={}", val_diff);
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
        methanol.set_2B_assoc_type(2899.5, 0.035176);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
        let mut methanol = PcSaftPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_3Bm_assoc_type(2899.5, 0.035176, 0.0);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
        let mut methanol = PcSaftPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_3B_assoc_type(2899.5, 0.035176);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
        let mut methanol = PcSaftPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_3Bm_assoc_type(2899.5, 0.035176, 1.0);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
    }
}
