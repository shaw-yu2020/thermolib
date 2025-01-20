use super::PcSaftErr;
use super::{A00, A01, A02, A03, A04, A05, A06, B00, B01, B02, B03, B04, B05, B06};
use super::{A10, A11, A12, A13, A14, A15, A16, B10, B11, B12, B13, B14, B15, B16};
use super::{A20, A21, A22, A23, A24, A25, A26, B20, B21, B22, B23, B24, B25, B26};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_6, PI};
/// PC-SAFT EOS
/// ```
/// use thermolib::PcSaftGlyPure;
/// let (m, sigma, epsilon) = (2.8611, 2.6826, 205.35); // SO2
/// let mut fluid = PcSaftGlyPure::new_fluid(m, sigma, epsilon);
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
/// let mut methanol = PcSaftGlyPure::new_fluid(1.5255, 3.23, 188.9);
/// methanol.set_2B_assoc_type(0.035176, 2899.5); // kappa_AB epsilon_AB
/// methanol.set_gly_params(1.0, 1.0, 1.0);
/// methanol.tp_flash(298.15, 0.1e6).unwrap();
/// assert_eq!(methanol.rho().unwrap().round(), 24676.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)]
pub struct PcSaftGlyPure {
    m: f64,
    epsilon: f64,
    sigma_plus: f64,
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
    // parameters for gly
    c0: f64,
    c1: f64,
    c2: f64,
}
impl PcSaftGlyPure {
    pub fn new_fluid(m: f64, sigma: f64, epsilon: f64) -> Self {
        let m1 = (m - 1.0) / m; // (m-1)/m
        let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
        Self {
            m,
            epsilon,
            sigma_plus: sigma.powi(3),
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
            // parameters for gly
            c0: 1.0,
            c1: -0.5,
            c2: 0.0,
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)]
impl PcSaftGlyPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_1_assoc_type(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type1 { X: 1.0 };
        self.kappa_AB_plus = kappa_AB * self.sigma_plus;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_2B_assoc_type(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type2B { X: 1.0 };
        self.kappa_AB_plus = kappa_AB * self.sigma_plus;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_3B_assoc_type(&mut self, kappa_AB: f64, epsilon_AB: f64) {
        self.assoc_type = AssocType::Type3B { XA: 1.0 };
        self.kappa_AB_plus = kappa_AB * self.sigma_plus;
        self.epsilon_AB = epsilon_AB;
    }
    pub fn set_aly_lee_cp0(&mut self, B: f64, C: f64, D: f64, E: f64, F: f64) {
        self.mB = B / R - 1.0; // mB = B/R -1
        self.mC = C / R; // mC = C/R
        self.mD = D; // mD = D
        self.mE = E / R; // mE = E/R
        self.mF = F; // mF = F
    }
    pub fn set_gly_params(&mut self, c0: f64, c1: f64, c2: f64) {
        self.c0 = c0;
        self.c1 = -2.0 * c0 + 1.5 * c1;
        self.c2 = c0 - 1.5 * c1 + 0.5 * c2;
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
            / ((FRAC_PI_6 * self.m * self.sigma_plus)
                * (1.0 - 0.12 * (-3.0 * self.epsilon / temp_c).exp()).powi(3));
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
        Err(anyhow!(PcSaftErr::NotConvForC))
    }
    pub fn t_flash(&mut self, temp: f64) -> anyhow::Result<()> {
        let d3 = self.sigma_plus * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp()).powi(3);
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
            return Err(anyhow!(PcSaftErr::NotConvForT));
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
            return Err(anyhow!(PcSaftErr::NotConvForT));
        }
        let pl_limit = self.calc_p(temp, rhol_num);
        if pl_limit > pv_limit {
            return Err(anyhow!(PcSaftErr::NotConvForT));
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
            Err(anyhow!(PcSaftErr::NotConvForT))
        } else {
            self.is_single_phase = false;
            self.temp = temp;
            self.rhov_num = rhov_num;
            self.rhol_num = rhol_num;
            Ok(())
        }
    }
    pub fn tp_flash(&mut self, temp: f64, p: f64) -> anyhow::Result<()> {
        let d3 = self.sigma_plus * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp()).powi(3);
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
    pub fn w(&mut self, M: f64) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok((self.calc_w2(self.temp, self.rho_num) / M).sqrt())
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
                let d3 = self.sigma_plus * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp()).powi(3);
                // Iteration from gas phase: eta = 1E-10
                self.calc_density(t, p, 1E-10 / (FRAC_PI_6 * self.m * d3))
            })
            .collect()
    }
    pub fn vec_tp_flash_l(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
        temp.into_iter()
            .zip(pres)
            .map(|(t, p)| {
                let d3 = self.sigma_plus * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp()).powi(3);
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
impl PcSaftGlyPure {
    fn calc_p(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp * rho_num * (1.0 + self.calc_rT0D1(temp, rho_num))
    }
    fn calc_Dp_Drho_T(&mut self, temp: f64, rho_num: f64) -> f64 {
        (FRAC_RE30_NA * temp)
            * (1.0 + 2.0 * self.calc_rT0D1(temp, rho_num) + self.calc_rT0D2(temp, rho_num))
    }
    fn calc_D2p_Drho2_T(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp / rho_num
            * (2.0 * self.calc_rT0D1(temp, rho_num)
                + 4.0 * self.calc_rT0D2(temp, rho_num)
                + self.calc_rT0D3(temp, rho_num))
    }
    fn calc_D3p_Drho3_T(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp / rho_num.powi(2)
            * (6.0 * self.calc_rT0D2(temp, rho_num)
                + 6.0 * self.calc_rT0D3(temp, rho_num)
                + self.calc_rT0D4(temp, rho_num))
    }
    fn calc_D2p_DTrho(&mut self, temp: f64, rho_num: f64) -> f64 {
        (FRAC_RE30_NA * temp)
            * (1.0
                + 2.0 * self.calc_rT0D1(temp, rho_num)
                + self.calc_rT0D2(temp, rho_num)
                + 2.0 * self.calc_rT1D1(temp, rho_num)
                + self.calc_rT1D2(temp, rho_num))
    }
    fn calc_D3p_DTrho2(&mut self, temp: f64, rho_num: f64) -> f64 {
        FRAC_RE30_NA * temp / rho_num
            * (2.0 * self.calc_rT0D1(temp, rho_num)
                + 4.0 * self.calc_rT0D2(temp, rho_num)
                + self.calc_rT0D3(temp, rho_num)
                + 2.0 * self.calc_rT1D1(temp, rho_num)
                + 4.0 * self.calc_rT1D2(temp, rho_num)
                + self.calc_rT1D3(temp, rho_num))
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
impl PcSaftGlyPure {
    fn set_temperature_and_number_density(&mut self, temp: f64, rho_num: f64) {
        if temp != self.temp {
            self.temp = temp;
            self.m2e1s3 = self.m.powi(2) * (self.epsilon / temp) * self.sigma_plus;
            self.m2e2s3 = self.m.powi(2) * (self.epsilon / temp).powi(2) * self.sigma_plus;
        } else if rho_num != self.rho_num {
        } else {
            return;
        }
        self.rho_num = rho_num;
        let d = 1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp();
        let d1 = -0.36 * (-3.0 * self.epsilon / temp).exp() * self.epsilon / temp;
        let d2 = d1 * (3.0 * self.epsilon / temp - 2.0);
        self.eta = FRAC_PI_6 * rho_num * self.m * self.sigma_plus * d.powi(3);
        self.eta1 = FRAC_PI_2 * rho_num * self.m * self.sigma_plus * d.powi(2) * d1;
        self.eta2 = PI * rho_num * self.m * self.sigma_plus * d * d1.powi(2)
            + FRAC_PI_2 * rho_num * self.m * self.sigma_plus * d.powi(2) * d2;
        if let AssocType::Type0 = self.assoc_type {
        } else {
            let t = self.tT0D0();
            self.assoc_type.set_t(t);
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
}
impl AssocType {
    fn set_t(&mut self, t: f64) {
        match self {
            AssocType::Type1 { X } | AssocType::Type2B { X } => {
                *X = (-1.0 + (1.0 + 4.0 * t).sqrt()) / (2.0 * t);
            }
            AssocType::Type3B { XA } => {
                *XA = (-(1.0 - t) + ((1.0 - t).powi(2) + 8.0 * t).sqrt()) / (4.0 * t);
            }

            _ => (),
        }
    }
    fn t1(&self) -> f64 {
        match self {
            AssocType::Type1 { X } | AssocType::Type2B { X } => X.powi(3) / (X - 2.0),
            AssocType::Type3B { XA } => {
                (XA * (2.0 * XA - 1.0)).powi(2) / (2.0 * XA.powi(2) - 4.0 * XA + 1.0)
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
            _ => 0.0,
        }
    }
}
#[allow(non_snake_case)]
impl PcSaftGlyPure {
    fn assocT0D0(&self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type1 { X } => X.ln() - X / 2.0 + 0.5,
            AssocType::Type2B { X } => 2.0 * X.ln() - X + 1.0,
            AssocType::Type3B { XA } => 2.0 * XA.ln() + (2.0 * XA - 1.0).ln() - 2.0 * XA + 2.0,
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
        }
    }
}
#[allow(non_snake_case)]
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
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
impl PcSaftGlyPure {
    fn DeltaT0D0(&mut self) -> f64 {
        self.giiT0D0_gly() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D1(&mut self) -> f64 {
        self.giiT0D1_gly() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D2(&mut self) -> f64 {
        self.giiT0D2_gly() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D3(&self) -> f64 {
        self.giiT0D3_gly() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT0D4(&self) -> f64 {
        self.giiT0D4_gly() * ((self.epsilon_AB / self.temp).exp() - 1.0)
    }
    fn DeltaT1D0(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D0_gly() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D0_gly() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT1D1(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D1_gly() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D1_gly() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT1D2(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D2_gly() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D2_gly() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT1D3(&self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT1D3_gly() * (epsilon_AB_T.exp() - 1.0)
            - self.giiT0D3_gly() * epsilon_AB_T.exp() * epsilon_AB_T
    }
    fn DeltaT2D0(&mut self) -> f64 {
        let epsilon_AB_T = self.epsilon_AB / self.temp;
        self.giiT2D0_gly() * (epsilon_AB_T.exp() - 1.0)
            - 2.0 * self.giiT1D0_gly() * epsilon_AB_T.exp() * epsilon_AB_T
            + self.giiT0D0_gly() * epsilon_AB_T.exp() * epsilon_AB_T * (epsilon_AB_T + 2.0)
    }
}
#[allow(non_snake_case)]
impl PcSaftGlyPure {
    fn giiT0D0_gly(&mut self) -> f64 {
        (1.0 - self.eta).powi(3).recip()
            * (self.c0 + self.c1 * self.eta + self.c2 * self.eta.powi(2))
    }
    fn giiT0D1_gly(&mut self) -> f64 {
        self.eta / (1.0 - self.eta).powi(4)
            * ((3.0 * self.c0 + self.c1)
                + (2.0 * self.c1 + 2.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn giiT0D2_gly(&mut self) -> f64 {
        self.eta.powi(2) * 2.0 / (1.0 - self.eta).powi(5)
            * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                + (3.0 * self.c1 + 4.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn giiT0D3_gly(&self) -> f64 {
        self.eta.powi(3) * 6.0 / (1.0 - self.eta).powi(6)
            * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                + (4.0 * self.c1 + 6.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn giiT0D4_gly(&self) -> f64 {
        self.eta.powi(4) * 24.0 / (1.0 - self.eta).powi(7)
            * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                + (5.0 * self.c1 + 8.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn giiT1D0_gly(&self) -> f64 {
        self.eta1 / (1.0 - self.eta).powi(4)
            * ((3.0 * self.c0 + self.c1)
                + (2.0 * self.c1 + 2.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn giiT1D1_gly(&self) -> f64 {
        self.eta1 / (1.0 - self.eta).powi(4)
            * ((3.0 * self.c0 + self.c1)
                + (2.0 * self.c1 + 2.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
            + self.eta1 * self.eta * 2.0 / (1.0 - self.eta).powi(5)
                * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                    + (3.0 * self.c1 + 4.0 * self.c2) * self.eta
                    + self.c2 * self.eta.powi(2))
    }
    fn giiT1D2_gly(&self) -> f64 {
        2.0 * self.eta1 * self.eta * 2.0 / (1.0 - self.eta).powi(5)
            * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                + (3.0 * self.c1 + 4.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
            + self.eta1 * self.eta.powi(2) * 6.0 / (1.0 - self.eta).powi(6)
                * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                    + (4.0 * self.c1 + 6.0 * self.c2) * self.eta
                    + self.c2 * self.eta.powi(2))
    }
    fn giiT1D3_gly(&self) -> f64 {
        3.0 * self.eta1 * self.eta.powi(2) * 6.0 / (1.0 - self.eta).powi(6)
            * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                + (4.0 * self.c1 + 6.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
            + self.eta1 * self.eta.powi(3) * 24.0 / (1.0 - self.eta).powi(7)
                * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                    + (5.0 * self.c1 + 8.0 * self.c2) * self.eta
                    + self.c2 * self.eta.powi(2))
    }
    fn giiT2D0_gly(&self) -> f64 {
        self.eta2 / (1.0 - self.eta).powi(4)
            * ((3.0 * self.c0 + self.c1)
                + (2.0 * self.c1 + 2.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
            + self.eta1.powi(2) * 2.0 / (1.0 - self.eta).powi(5)
                * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                    + (3.0 * self.c1 + 4.0 * self.c2) * self.eta
                    + self.c2 * self.eta.powi(2))
    }
}
const FRAC_RE30_NA: f64 = R / FRAC_NA_1E30;
const FRAC_NA_1E30: f64 = 6.02214076E-7; // const NA: f64 = 6.02214076E23;
const R: f64 = 8.314462618;
impl PcSaftGlyPure {
    pub fn check_derivatives(&mut self, print_val: bool) {
        let (t, d) = (self.temp, self.rho_num);
        if print_val {
            println!("[rT0D0 == rT0D0] calc_rT0D0() ={}", self.calc_rT0D0(t, d));
        }
        let compare_val = |val_calc: f64, val_diff: f64| {
            assert_eq!(
                &val_calc.abs().to_string()[0..11],
                &val_diff.abs().to_string()[0..11]
            )
        };
        // derivative for density
        let val_calc = self.calc_rT0D1(t, d) / d;
        let val_diff = romberg_diff(|dx: f64| self.calc_rT0D0(t, dx), d);
        if print_val {
            println!("[rT0D1 == rT0D1] calc_rT0D1() ={}", val_calc);
            println!("[rT0D0 -> rT0D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density
        let val_calc = self.calc_rT0D2(t, d) / d.powi(2);
        let val_diff = romberg_diff(|dx: f64| self.calc_rT0D1(t, dx) / dx, d);
        if print_val {
            println!("[rT0D2 == rT0D2] calc_rT0D2() ={}", val_calc);
            println!("[rT0D1 -> rT0D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density
        let val_calc = self.calc_rT0D3(t, d) / d.powi(3);
        let val_diff = romberg_diff(|dx: f64| self.calc_rT0D2(t, dx) / dx.powi(2), d);
        if print_val {
            println!("[rT0D3 == rT0D3] calc_rT0D3() ={}", val_calc);
            println!("[rT0D2 -> rT0D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density+density
        let val_calc = self.calc_rT0D4(t, d) / d.powi(4);
        let val_diff = romberg_diff(|dx: f64| self.calc_rT0D3(t, dx) / dx.powi(3), d);
        if print_val {
            println!("[rT0D4 == rT0D4] calc_rT0D4() ={}", val_calc);
            println!("[rT0D3 -> rT0D4] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature
        let val_calc = self.calc_rT1D0(t, d) / t;
        let val_diff = romberg_diff(|tx: f64| self.calc_rT0D0(tx, d), t);
        if print_val {
            println!("[rT1D0 == rT1D0] calc_rT1D0() ={}", val_calc);
            println!("[rT0D0 -> rT1D0] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density
        let val_calc = self.calc_rT1D1(t, d) / t / d;
        let val_diff = romberg_diff(|dx: f64| self.calc_rT1D0(t, dx) / t, d);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_rT1D1() ={}", val_calc);
            println!("[rT1D0 -> rT1D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.calc_rT0D1(tx, d) / d, t);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_rT1D1() ={}", val_calc);
            println!("[rT0D1 -> rT1D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density
        let val_calc = self.calc_rT1D2(t, d) / t / d.powi(2);
        let val_diff = romberg_diff(|dx: f64| self.calc_rT1D1(t, dx) / t / dx, d);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_rT1D2() ={}", val_calc);
            println!("[rT1D1 -> rT1D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.calc_rT0D2(tx, d) / d.powi(2), t);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_rT1D2() ={}", val_calc);
            println!("[rT0D2 -> rT1D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density+density
        let val_calc = self.calc_rT1D3(t, d) / t / d.powi(3);
        let val_diff = romberg_diff(|dx: f64| self.calc_rT1D2(t, dx) / t / dx.powi(2), d);
        if print_val {
            println!("[rT1D3 == rT1D3] calc_rT1D3() ={}", val_calc);
            println!("[rT1D2 -> rT1D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.calc_rT0D3(tx, d) / d.powi(3), t);
        if print_val {
            println!("[rT1D3 == rT1D3] calc_rT1D3() ={}", val_calc);
            println!("[rT0D3 -> rT1D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+temperature
        let val_calc = self.calc_rT2D0(t, d) / t.powi(2);
        let val_diff = romberg_diff(|tx: f64| self.calc_rT1D0(tx, d) / tx, t);
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
        let mut fluid = PcSaftGlyPure::new_fluid(m, sigma, epsilon);
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
        let mut methanol = PcSaftGlyPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_2B_assoc_type(0.035176, 2899.5);
        methanol.set_gly_params(1.0, 1.0, 1.0);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
        let mut methanol = PcSaftGlyPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_3B_assoc_type(0.035176, 2899.5);
        methanol.set_gly_params(1.0, 1.0, 1.0);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
    }
}
