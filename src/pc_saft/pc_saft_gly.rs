use super::{CTerm, I1Term, I2Term, PcSaftErr};
use crate::algorithms::{brent_zero, romberg_diff};
use crate::f64consts::{FRAC_NA_1E30, FRAC_PI_2, FRAC_PI_6, PI, R};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
/// PC-SAFT EOS
/// ```
/// use thermolib::PcSaftGlyPure;
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
pub struct PcSaftGlyPure {
    m: f64,
    sigma_plus: f64,
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
    kappa_ab_plus: f64,
    epsilon_ab: f64,
    xa: f64,
    // modified aly_lee_cv0
    cv_b: f64, // b=B/R-1
    cv_c: f64, // c=C/R
    cv_d: f64, // d=D
    cv_e: f64, // e=E/R
    cv_f: f64, // f=F
    // cached variables
    gii_t0d0: (f64, f64),
    gii_t0d1: (f64, f64),
    gii_t0d2: (f64, f64),
    i1: I1Term, // I1Term
    i2: I2Term, // I2Term
    c: CTerm,   // Cterm
    // parameters for gly
    c0: f64,
    c1: f64,
    c2: f64,
}
impl PcSaftGlyPure {
    pub fn new_fluid(m: f64, sigma: f64, epsilon: f64) -> Self {
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
            kappa_ab_plus: 0.0,
            epsilon_ab: 0.0,
            xa: 1.0,
            // modified aly_lee_cv0
            cv_b: 1.5,
            cv_c: 0.0,
            cv_d: 1.0,
            cv_e: 0.0,
            cv_f: 1.0,
            // cached variables
            gii_t0d0: (0.0, 0.0),
            gii_t0d1: (0.0, 0.0),
            gii_t0d2: (0.0, 0.0),
            i1: I1Term::new(m), // I1Term
            i2: I2Term::new(m), // I1Term
            c: CTerm::new(m),   // CTerm
            // parameters for gly
            c0: 1.0,
            c1: -0.5,
            c2: 0.0,
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
impl PcSaftGlyPure {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self::new_fluid(m, sigma, epsilon)
    }
    pub fn set_2b_assoc_type(&mut self, kappa_ab: f64, epsilon_ab: f64) {
        self.assoc_type = AssocType::Type2B;
        self.kappa_ab_plus = kappa_ab * self.sigma_plus;
        self.epsilon_ab = epsilon_ab;
    }
    pub fn set_3b_assoc_type(&mut self, kappa_ab: f64, epsilon_ab: f64) {
        self.assoc_type = AssocType::Type3B;
        self.kappa_ab_plus = kappa_ab * self.sigma_plus;
        self.epsilon_ab = epsilon_ab;
    }
    pub fn set_aly_lee_cp0(&mut self, b: f64, c: f64, d: f64, e: f64, f: f64) {
        self.cv_b = b / R - 1.0; // B/R -1
        self.cv_c = c / R; // C/R
        self.cv_d = d; // D
        self.cv_e = e / R; // E/R
        self.cv_f = f; // F
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
        let mut temp_c: f64 = 1000.0;
        let mut dens_c: f64 = 1E-10
            / ((FRAC_PI_6 * self.m * self.sigma_plus)
                * (1.0 - 0.12 * (-3.0 * self.epsilon / temp_c).exp()).powi(3));
        // Define variables
        let mut drho = self.calc_p_drho(temp_c, dens_c);
        let mut drho2 = self.calc_p_drho2(temp_c, dens_c);
        let (mut drho3, mut dt_drho, mut dt_drho2);
        for _i in 1..100000 {
            drho3 = self.calc_p_drho3(temp_c, dens_c);
            dt_drho = self.calc_p_dt_drho(temp_c, dens_c);
            dt_drho2 = self.calc_p_dt_drho2(temp_c, dens_c);
            temp_c -= (drho * drho3 - drho2 * drho2) / (dt_drho * drho3 - dt_drho2 * drho2);
            dens_c -= (drho * dt_drho2 - drho2 * dt_drho) / (drho2 * dt_drho2 - drho3 * dt_drho);
            drho = self.calc_p_drho(temp_c, dens_c);
            drho2 = self.calc_p_drho2(temp_c, dens_c);
            if drho.abs() < 1E3 && drho2.abs() < 1E6 {
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
        let mut drhov = self.calc_p_drho(temp, rhov_num);
        for _i in 1..100000 {
            if drhov.abs() < 10.0 {
                break;
            } else {
                rhov_num -= drhov / self.calc_p_drho2(temp, rhov_num);
                drhov = self.calc_p_drho(temp, rhov_num);
            }
        }
        let pv_limit = self.calc_p(temp, rhov_num);
        if pv_limit.is_sign_negative() {
            return Err(anyhow!(PcSaftErr::NotConvForT));
        }
        // Liquid phase: eta = 0.5
        let rhol_num_guess = 0.5 / (FRAC_PI_6 * self.m * d3);
        let mut rhol_num = rhol_num_guess;
        let mut drhol = self.calc_p_drho(temp, rhol_num);
        for _i in 1..100000 {
            if drhol.abs() < 1.0 {
                break;
            } else {
                rhol_num -= drhol / self.calc_p_drho2(temp, rhol_num);
                drhol = self.calc_p_drho(temp, rhol_num);
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
    pub fn t(&self) -> anyhow::Result<f64> {
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
    pub fn t_s(&self) -> anyhow::Result<f64> {
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
    pub fn vec_w(&mut self, temp: Vec<f64>, dens_num: Vec<f64>, molar_mass: f64) -> Vec<f64> {
        temp.into_iter()
            .zip(dens_num)
            .map(|(t, d)| (self.calc_w2(t, d) / molar_mass).sqrt()) // m/s
            .collect()
    }
}
impl PcSaftGlyPure {
    fn calc_p(&mut self, temp: f64, rho_num: f64) -> f64 {
        R / FRAC_NA_1E30 * temp * rho_num * (1.0 + self.r_t0d1(temp, rho_num))
    }
    fn calc_p_drho(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (1.0 + 2.0 * self.r_t0d1(temp, rho_num) + self.r_t0d2(temp, rho_num))
    }
    fn calc_p_drho2(&mut self, temp: f64, rho_num: f64) -> f64 {
        R / FRAC_NA_1E30 * temp / rho_num
            * (2.0 * self.r_t0d1(temp, rho_num)
                + 4.0 * self.r_t0d2(temp, rho_num)
                + self.r_t0d3(temp, rho_num))
    }
    fn calc_p_drho3(&mut self, temp: f64, rho_num: f64) -> f64 {
        R / FRAC_NA_1E30 * temp / rho_num.powi(2)
            * (6.0 * self.r_t0d2(temp, rho_num)
                + 6.0 * self.r_t0d3(temp, rho_num)
                + self.r_t0d4(temp, rho_num))
    }
    fn calc_p_dt_drho(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (1.0
                + 2.0 * self.r_t0d1(temp, rho_num)
                + self.r_t0d2(temp, rho_num)
                + 2.0 * self.r_t1d1(temp, rho_num)
                + self.r_t1d2(temp, rho_num))
    }
    fn calc_p_dt_drho2(&mut self, temp: f64, rho_num: f64) -> f64 {
        R / FRAC_NA_1E30 * temp / rho_num
            * (2.0 * self.r_t0d1(temp, rho_num)
                + 4.0 * self.r_t0d2(temp, rho_num)
                + self.r_t0d3(temp, rho_num)
                + 2.0 * self.r_t1d1(temp, rho_num)
                + 4.0 * self.r_t1d2(temp, rho_num)
                + self.r_t1d3(temp, rho_num))
    }
    fn calc_w2(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * temp
            * ((1.0 + 2.0 * self.r_t0d1(temp, rho_num) + self.r_t0d2(temp, rho_num))
                + (1.0 + self.r_t0d1(temp, rho_num) + self.r_t1d1(temp, rho_num)).powi(2)
                    / (self.calc_ideal_cv(temp)
                        - 2.0 * self.r_t1d0(temp, rho_num)
                        - self.r_t2d0(temp, rho_num)))
    }
    fn calc_cv(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * (self.calc_ideal_cv(temp)
            - 2.0 * self.r_t1d0(temp, rho_num)
            - self.r_t2d0(temp, rho_num))
    }
    fn calc_cp(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * (self.calc_ideal_cv(temp)
            - 2.0 * self.r_t1d0(temp, rho_num)
            - self.r_t2d0(temp, rho_num)
            + (1.0 + self.r_t0d1(temp, rho_num) + self.r_t1d1(temp, rho_num)).powi(2)
                / (1.0 + 2.0 * self.r_t0d1(temp, rho_num) + self.r_t0d2(temp, rho_num)))
    }
    fn calc_cp_res(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * ((-1.0 - 2.0 * self.r_t1d0(temp, rho_num) - self.r_t2d0(temp, rho_num))
            + (1.0 + self.r_t0d1(temp, rho_num) + self.r_t1d1(temp, rho_num)).powi(2)
                / (1.0 + 2.0 * self.r_t0d1(temp, rho_num) + self.r_t0d2(temp, rho_num)))
    }
    fn calc_h_res(&mut self, temp: f64, rho_num: f64) -> f64 {
        R * temp * (-self.r_t1d0(temp, rho_num) + self.r_t0d1(temp, rho_num))
    }
    fn calc_lnphi(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.r_t0d0(temp, rho_num) + self.r_t0d1(temp, rho_num)
            - (1.0 + self.r_t0d1(temp, rho_num)).ln()
    }
    fn calc_density(&mut self, temp: f64, p: f64, rho_num_guess: f64) -> f64 {
        let mut rho_num = rho_num_guess;
        let (mut p_diff, mut val_p_drho, mut rho_num_diff);
        for _i in 1..10000 {
            p_diff = self.calc_p(temp, rho_num) - p;
            if p_diff.abs() < f64::EPSILON {
                return rho_num;
            }
            val_p_drho = self.calc_p_drho(temp, rho_num);
            if val_p_drho.is_sign_negative() {
                return f64::NAN;
            }
            rho_num_diff = p_diff / val_p_drho;
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
            let t = self.t_t0d0();
            match self.assoc_type {
                AssocType::Type2B => {
                    self.xa = (-1.0 + (1.0 + 4.0 * t).sqrt()) / (2.0 * t);
                }
                AssocType::Type3B => {
                    self.xa = (-(1.0 - t) + ((1.0 - t).powi(2) + 8.0 * t).sqrt()) / (4.0 * t);
                }
                _ => (),
            }
        }
    }
    fn calc_ideal_cv(&mut self, temp: f64) -> f64 {
        self.cv_b
            + self.cv_c * (self.cv_d / temp / (self.cv_d / temp).sinh()).powi(2)
            + self.cv_e * (self.cv_f / temp / (self.cv_f / temp).cosh()).powi(2)
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d0() + self.disp_t0d0() + self.assoc_t0d0()
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d1() + self.disp_t0d1() + self.assoc_t0d1()
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d2() + self.disp_t0d2() + self.assoc_t0d2()
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d3() + self.disp_t0d3() + self.assoc_t0d3()
    }
    fn r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t0d4() + self.disp_t0d4() + self.assoc_t0d4()
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d0() + self.disp_t1d0() + self.assoc_t1d0()
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d1() + self.disp_t1d1() + self.assoc_t1d1()
    }
    fn r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d2() + self.disp_t1d2() + self.assoc_t1d2()
    }
    fn r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t1d3() + self.disp_t1d3() + self.assoc_t1d3()
    }
    fn r_t2d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(temp, rho_num);
        self.hc_t2d0() + self.disp_t2d0() + self.assoc_t2d0()
    }
}
impl PcSaftGlyPure {
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
        let gii_t1d0: f64 = self.gii_t1d0();
        self.m * self.hs_t1d2()
            + (1.0 - self.m)
                * (2.0 / self.gii_t0d0().powi(3) * gii_t1d0 * self.gii_t0d1().powi(2)
                    - 2.0 / self.gii_t0d0().powi(2) * self.gii_t1d1() * self.gii_t0d1()
                    - gii_t1d0 * self.gii_t0d2() / self.gii_t0d0().powi(2)
                    + self.gii_t1d2() / self.gii_t0d0())
    }
    fn hc_t1d3(&mut self) -> f64 {
        let gii_t1d0: f64 = self.gii_t1d0();
        let gii_t1d1: f64 = self.gii_t1d1();
        self.m * self.hs_t1d3()
            + (1.0 - self.m)
                * (-6.0 / self.gii_t0d0().powi(4) * gii_t1d0 * self.gii_t0d1().powi(3)
                    + 6.0 / self.gii_t0d0().powi(3) * gii_t1d1 * self.gii_t0d1().powi(2)
                    + 6.0 / self.gii_t0d0().powi(3) * gii_t1d0 * self.gii_t0d1() * self.gii_t0d2()
                    - 3.0 / self.gii_t0d0().powi(2) * gii_t1d1 * self.gii_t0d2()
                    - 3.0 / self.gii_t0d0().powi(2) * self.gii_t1d2() * self.gii_t0d1()
                    - gii_t1d0 * self.gii_t0d3() / self.gii_t0d0().powi(2)
                    + self.gii_t1d3() / self.gii_t0d0())
    }
    fn hc_t2d0(&mut self) -> f64 {
        self.m * self.hs_t2d0()
            + (1.0 - self.m)
                * (-self.gii_t1d0().powi(2) / self.gii_t0d0().powi(2)
                    + self.gii_t2d0() / self.gii_t0d0())
    }
}
impl PcSaftGlyPure {
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
        self.eta1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
    }
    fn hs_t1d1(&self) -> f64 {
        self.eta1 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
    fn hs_t1d2(&self) -> f64 {
        2.0 * self.eta1 * self.eta * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
            + self.eta1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn hs_t1d3(&self) -> f64 {
        3.0 * self.eta1 * self.eta.powi(2) * (36.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(5)
            + self.eta1 * self.eta.powi(3) * (168.0 - 48.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn hs_t2d0(&self) -> f64 {
        self.eta2 * (4.0 - 2.0 * self.eta) / (1.0 - self.eta).powi(3)
            + self.eta1.powi(2) * (10.0 - 4.0 * self.eta) / (1.0 - self.eta).powi(4)
    }
}
impl PcSaftGlyPure {
    fn gii_t0d0(&mut self) -> f64 {
        if self.eta != self.gii_t0d0.0 {
            self.gii_t0d0 = (self.eta, (1.0 - 0.5 * self.eta) / (1.0 - self.eta).powi(3))
        }
        self.gii_t0d0.1
    }
    fn gii_t0d1(&mut self) -> f64 {
        if self.eta != self.gii_t0d1.0 {
            self.gii_t0d1 = (
                self.eta,
                self.eta * (2.5 - self.eta) / (1.0 - self.eta).powi(4),
            )
        }
        self.gii_t0d1.1
    }
    fn gii_t0d2(&mut self) -> f64 {
        if self.eta != self.gii_t0d2.0 {
            self.gii_t0d2 = (
                self.eta,
                self.eta.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5),
            )
        }
        self.gii_t0d2.1
    }
    fn gii_t0d3(&self) -> f64 {
        self.eta.powi(3) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn gii_t0d4(&self) -> f64 {
        self.eta.powi(4) * (240.0 - 60.0 * self.eta) / (1.0 - self.eta).powi(7)
    }
    fn gii_t1d0(&self) -> f64 {
        self.eta1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
    }
    fn gii_t1d1(&self) -> f64 {
        self.eta1 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
            + self.eta1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
    fn gii_t1d2(&self) -> f64 {
        2.0 * self.eta1 * self.eta * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
            + self.eta1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
    }
    fn gii_t1d3(&self) -> f64 {
        3.0 * self.eta1 * self.eta.powi(2) * (42.0 - 12.0 * self.eta) / (1.0 - self.eta).powi(6)
            + self.eta1 * self.eta.powi(3) * (240.0 - 60.0 * self.eta) / (1.0 - self.eta).powi(7)
    }
    fn gii_t2d0(&self) -> f64 {
        self.eta2 * (2.5 - self.eta) / (1.0 - self.eta).powi(4)
            + self.eta1.powi(2) * (9.0 - 3.0 * self.eta) / (1.0 - self.eta).powi(5)
    }
}
impl PcSaftGlyPure {
    fn disp_t0d0(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * self.i1t0d0(self.eta)
                + self.m * self.m2e2s3 * self.c1_t0d0() * self.i2t0d0())
    }
    fn disp_t0d1(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (self.i1t0d0(self.eta) + self.i1t0d1(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1_t0d0() * self.i2t0d0()
                        + self.c1_t0d1() * self.i2t0d0()
                        + self.c1_t0d0() * self.i2t0d1()))
    }
    fn disp_t0d2(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (2.0 * self.i1t0d1(self.eta) + self.i1t0d2(self.eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * self.c1_t0d1() * self.i2t0d0()
                        + 2.0 * self.c1_t0d0() * self.i2t0d1()
                        + self.c1_t0d2() * self.i2t0d0()
                        + 2.0 * self.c1_t0d1() * self.i2t0d1()
                        + self.c1_t0d0() * self.i2t0d2()))
    }
    fn disp_t0d3(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (3.0 * self.i1t0d2(self.eta) + self.i1t0d3(self.eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * self.c1_t0d2() * self.i2t0d0()
                        + 6.0 * self.c1_t0d1() * self.i2t0d1()
                        + 3.0 * self.c1_t0d0() * self.i2t0d2()
                        + self.c1_t0d3(self.eta) * self.i2t0d0()
                        + 3.0 * self.c1_t0d2() * self.i2t0d1()
                        + 3.0 * self.c1_t0d1() * self.i2t0d2()
                        + self.c1_t0d0() * self.i2t0d3(self.eta)))
    }
    fn disp_t0d4(&mut self) -> f64 {
        let c1_t0d3 = self.c1_t0d3(self.eta);
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (4.0 * self.i1t0d3(self.eta) + self.i1t0d4(self.eta))
                + self.m
                    * self.m2e2s3
                    * (4.0 * c1_t0d3 * self.i2t0d0()
                        + 12.0 * self.c1_t0d2() * self.i2t0d1()
                        + 12.0 * self.c1_t0d1() * self.i2t0d2()
                        + 4.0 * self.c1_t0d0() * self.i2t0d3(self.eta)
                        + self.c1_t0d4(self.eta) * self.i2t0d0()
                        + 4.0 * c1_t0d3 * self.i2t0d1()
                        + 6.0 * self.c1_t0d2() * self.i2t0d2()
                        + 4.0 * self.c1_t0d1() * self.i2t0d3(self.eta)
                        + self.c1_t0d0() * self.i2t0d4(self.eta)))
    }
    fn disp_t1d0(&mut self) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (self.i1t1d0(self.eta) - self.i1t0d0(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1_t1d0(self.eta) * self.i2t0d0()
                        + self.c1_t0d0() * self.i2t1d0(self.eta)
                        - 2.0 * self.c1_t0d0() * self.i2t0d0()))
    }
    fn disp_t1d1(&mut self) -> f64 {
        let c1_t1d0 = self.c1_t1d0(self.eta);
        let i2_t1d0 = self.i2t1d0(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (self.i1t1d0(self.eta) - self.i1t0d0(self.eta) + self.i1t1d1(self.eta)
                    - self.i1t0d1(self.eta))
                + self.m
                    * self.m2e2s3
                    * (c1_t1d0 * self.i2t0d0() + self.c1_t0d0() * i2_t1d0
                        - 2.0 * self.c1_t0d0() * self.i2t0d0()
                        + self.c1_t1d1(self.eta) * self.i2t0d0()
                        + self.c1_t0d1() * i2_t1d0
                        - 2.0 * self.c1_t0d1() * self.i2t0d0()
                        + c1_t1d0 * self.i2t0d1()
                        + self.c1_t0d0() * self.i2t1d1(self.eta)
                        - 2.0 * self.c1_t0d0() * self.i2t0d1()))
    }
    fn disp_t1d2(&mut self) -> f64 {
        let c1_t1d0 = self.c1_t1d0(self.eta);
        let c1_t1d1 = self.c1_t1d1(self.eta);
        let i2_t1d0 = self.i2t1d0(self.eta);
        let i2_t1d1 = self.i2t1d1(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (2.0 * self.i1t1d1(self.eta) - 2.0 * self.i1t0d1(self.eta)
                    + self.i1t1d2(self.eta)
                    - self.i1t0d2(self.eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * c1_t1d1 * self.i2t0d0() + 2.0 * self.c1_t0d1() * i2_t1d0
                        - 4.0 * self.c1_t0d1() * self.i2t0d0()
                        + 2.0 * c1_t1d0 * self.i2t0d1()
                        + 2.0 * self.c1_t0d0() * i2_t1d1
                        - 4.0 * self.c1_t0d0() * self.i2t0d1()
                        + self.c1_t1d2(self.eta) * self.i2t0d0()
                        + self.c1_t0d2() * i2_t1d0
                        - 2.0 * self.c1_t0d2() * self.i2t0d0()
                        + 2.0 * c1_t1d1 * self.i2t0d1()
                        + 2.0 * self.c1_t0d1() * i2_t1d1
                        - 4.0 * self.c1_t0d1() * self.i2t0d1()
                        + c1_t1d0 * self.i2t0d2()
                        + self.c1_t0d0() * self.i2t1d2(self.eta)
                        - 2.0 * self.c1_t0d0() * self.i2t0d2()))
    }
    fn disp_t1d3(&mut self) -> f64 {
        let c1_t0d3 = self.c1_t0d3(self.eta);
        let c1_t1d0 = self.c1_t1d0(self.eta);
        let c1_t1d1 = self.c1_t1d1(self.eta);
        let c1_t1d2 = self.c1_t1d2(self.eta);
        let i2_t1d0 = self.i2t1d0(self.eta);
        let i2_t1d1 = self.i2t1d1(self.eta);
        let i2_t1d2 = self.i2t1d2(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (3.0 * self.i1t1d2(self.eta) - 3.0 * self.i1t0d2(self.eta)
                    + self.i1t1d3(self.eta)
                    - self.i1t0d3(self.eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * self.c1_t1d2(self.eta) * self.i2t0d0()
                        + 3.0 * self.c1_t0d2() * i2_t1d0
                        - 6.0 * self.c1_t0d2() * self.i2t0d0()
                        + 6.0 * c1_t1d1 * self.i2t0d1()
                        + 6.0 * self.c1_t0d1() * i2_t1d1
                        - 12.0 * self.c1_t0d1() * self.i2t0d1()
                        + 3.0 * c1_t1d0 * self.i2t0d2()
                        + 3.0 * self.c1_t0d0() * self.i2t1d2(self.eta)
                        - 6.0 * self.c1_t0d0() * self.i2t0d2()
                        + self.c1_t1d3(self.eta) * self.i2t0d0()
                        + c1_t0d3 * i2_t1d0
                        - 2.0 * c1_t0d3 * self.i2t0d0()
                        + 3.0 * c1_t1d2 * self.i2t0d1()
                        + 3.0 * self.c1_t0d2() * i2_t1d1
                        - 6.0 * self.c1_t0d2() * self.i2t0d1()
                        + 3.0 * c1_t1d1 * self.i2t0d2()
                        + 3.0 * self.c1_t0d1() * i2_t1d2
                        - 6.0 * self.c1_t0d1() * self.i2t0d2()
                        + c1_t1d0 * self.i2t0d3(self.eta)
                        + self.c1_t0d0() * self.i2t1d3(self.eta)
                        - 2.0 * self.c1_t0d0() * self.i2t0d3(self.eta)))
    }
    fn disp_t2d0(&mut self) -> f64 {
        let c1_t1d0 = self.c1_t1d0(self.eta);
        let i2_t1d0 = self.i2t1d0(self.eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (self.i1t2d0(self.eta) + 2.0 * self.i1t0d0(self.eta)
                    - 2.0 * self.i1t1d0(self.eta))
                + self.m
                    * self.m2e2s3
                    * (self.c1_t2d0(self.eta) * self.i2t0d0()
                        + self.c1_t0d0() * self.i2t2d0(self.eta)
                        + 6.0 * self.c1_t0d0() * self.i2t0d0()
                        + 2.0 * c1_t1d0 * i2_t1d0
                        - 4.0 * c1_t1d0 * self.i2t0d0()
                        - 4.0 * self.c1_t0d0() * i2_t1d0))
    }
}
impl PcSaftGlyPure {
    fn c1_t0d0(&mut self) -> f64 {
        self.c_t0d0().recip()
    }
    fn c1_t0d1(&mut self) -> f64 {
        -self.c_t0d1() / self.c_t0d0().powi(2)
    }
    fn c1_t0d2(&mut self) -> f64 {
        2.0 * self.c_t0d1().powi(2) / self.c_t0d0().powi(3) - self.c_t0d2() / self.c_t0d0().powi(2)
    }
    fn c1_t0d3(&mut self, eta: f64) -> f64 {
        -6.0 * self.c_t0d1().powi(3) / self.c_t0d0().powi(4)
            + 6.0 * self.c_t0d1() * self.c_t0d2() / self.c_t0d0().powi(3)
            - self.c_t0d3(eta) / self.c_t0d0().powi(2)
    }
    fn c1_t0d4(&mut self, eta: f64) -> f64 {
        24.0 * self.c_t0d1().powi(4) / self.c_t0d0().powi(5)
            - 36.0 * self.c_t0d1().powi(2) * self.c_t0d2() / self.c_t0d0().powi(4)
            + 8.0 * self.c_t0d1() * self.c_t0d3(eta) / self.c_t0d0().powi(3)
            + 6.0 * self.c_t0d2().powi(2) / self.c_t0d0().powi(3)
            - self.c_t0d4(eta) / self.c_t0d0().powi(2)
    }
    fn c1_t1d0(&mut self, eta: f64) -> f64 {
        -self.c_t1d0(eta) / self.c_t0d0().powi(2)
    }
    fn c1_t1d1(&mut self, eta: f64) -> f64 {
        2.0 * self.c_t1d0(eta) * self.c_t0d1() / self.c_t0d0().powi(3)
            - self.c_t1d1(eta) / self.c_t0d0().powi(2)
    }
    fn c1_t1d2(&mut self, eta: f64) -> f64 {
        let c_t1d0 = self.c_t1d0(eta);
        -6.0 * c_t1d0 * self.c_t0d1().powi(2) / self.c_t0d0().powi(4)
            + 4.0 * self.c_t1d1(eta) * self.c_t0d1() / self.c_t0d0().powi(3)
            + 2.0 * c_t1d0 * self.c_t0d2() / self.c_t0d0().powi(3)
            - self.c_t1d2(eta) / self.c_t0d0().powi(2)
    }
    fn c1_t1d3(&mut self, eta: f64) -> f64 {
        let c_t1d0 = self.c_t1d0(eta);
        let c_t1d1 = self.c_t1d1(eta);
        24.0 * c_t1d0 * self.c_t0d1().powi(3) / self.c_t0d0().powi(5)
            - 18.0 * c_t1d1 * self.c_t0d1().powi(2) / self.c_t0d0().powi(4)
            - 18.0 * c_t1d0 * self.c_t0d1() * self.c_t0d2() / self.c_t0d0().powi(4)
            + 6.0 * c_t1d1 * self.c_t0d2() / self.c_t0d0().powi(3)
            + 6.0 * self.c_t1d2(eta) * self.c_t0d1() / self.c_t0d0().powi(3)
            + 2.0 * c_t1d0 * self.c_t0d3(eta) / self.c_t0d0().powi(3)
            - self.c_t1d3(eta) / self.c_t0d0().powi(2)
    }
    fn c1_t2d0(&mut self, eta: f64) -> f64 {
        2.0 * self.c_t1d0(eta).powi(2) / self.c_t0d0().powi(3)
            - self.c_t2d0(eta) / self.c_t0d0().powi(2)
    }
}
impl PcSaftGlyPure {
    fn c_t0d0(&mut self) -> f64 {
        self.c.eta0(self.eta)
    }
    fn c_t0d1(&mut self) -> f64 {
        self.eta * self.c.eta1(self.eta)
    }
    fn c_t0d2(&mut self) -> f64 {
        self.eta.powi(2) * self.c.eta2(self.eta)
    }
    fn c_t0d3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.c.eta3(eta)
    }
    fn c_t0d4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.c.eta4(eta)
    }
    fn c_t1d0(&mut self, eta: f64) -> f64 {
        self.eta1 * self.c.eta1(eta)
    }
    fn c_t1d1(&mut self, eta: f64) -> f64 {
        self.eta1 * (self.c.eta1(eta) + self.eta * self.c.eta2(eta))
    }
    fn c_t1d2(&mut self, eta: f64) -> f64 {
        self.eta1 * self.eta * (2.0 * self.c.eta2(eta) + self.eta * self.c.eta3(eta))
    }
    fn c_t1d3(&mut self, eta: f64) -> f64 {
        self.eta1 * self.eta.powi(2) * (3.0 * self.c.eta3(eta) + self.eta * self.c.eta4(eta))
    }
    fn c_t2d0(&mut self, eta: f64) -> f64 {
        self.eta2 * self.c.eta1(eta) + self.eta1.powi(2) * self.c.eta2(eta)
    }
}
impl PcSaftGlyPure {
    fn i1t0d0(&mut self, eta: f64) -> f64 {
        self.i1.eta0(eta)
    }
    fn i1t0d1(&mut self, eta: f64) -> f64 {
        eta * self.i1.eta1(eta)
    }
    fn i1t0d2(&mut self, eta: f64) -> f64 {
        eta.powi(2) * self.i1.eta2(eta)
    }
    fn i1t0d3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.i1.eta3(eta)
    }
    fn i1t0d4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.i1.eta4(eta)
    }
    fn i1t1d0(&mut self, eta: f64) -> f64 {
        self.eta1 * self.i1.eta1(eta)
    }
    fn i1t1d1(&mut self, eta: f64) -> f64 {
        self.eta1 * (self.i1.eta1(eta) + self.eta * self.i1.eta2(eta))
    }
    fn i1t1d2(&mut self, eta: f64) -> f64 {
        self.eta1 * self.eta * (2.0 * self.i1.eta2(eta) + self.eta * self.i1.eta3(eta))
    }
    fn i1t1d3(&mut self, eta: f64) -> f64 {
        self.eta1 * self.eta.powi(2) * (3.0 * self.i1.eta3(eta) + self.eta * self.i1.eta4(eta))
    }
    fn i1t2d0(&mut self, eta: f64) -> f64 {
        self.eta2 * self.i1.eta1(eta) + self.eta1.powi(2) * self.i1.eta2(eta)
    }
}
impl PcSaftGlyPure {
    fn i2t0d0(&mut self) -> f64 {
        self.i2.eta0(self.eta)
    }
    fn i2t0d1(&mut self) -> f64 {
        self.eta * self.i2.eta1(self.eta)
    }
    fn i2t0d2(&mut self) -> f64 {
        self.eta.powi(2) * self.i2.eta2(self.eta)
    }
    fn i2t0d3(&mut self, eta: f64) -> f64 {
        eta.powi(3) * self.i2.eta3(eta)
    }
    fn i2t0d4(&mut self, eta: f64) -> f64 {
        eta.powi(4) * self.i2.eta4(eta)
    }
    fn i2t1d0(&mut self, eta: f64) -> f64 {
        self.eta1 * self.i2.eta1(eta)
    }
    fn i2t1d1(&mut self, eta: f64) -> f64 {
        self.eta1 * (self.i2.eta1(eta) + self.eta * self.i2.eta2(eta))
    }
    fn i2t1d2(&mut self, eta: f64) -> f64 {
        self.eta1 * self.eta * (2.0 * self.i2.eta2(eta) + self.eta * self.i2.eta3(eta))
    }
    fn i2t1d3(&mut self, eta: f64) -> f64 {
        self.eta1 * self.eta.powi(2) * (3.0 * self.i2.eta3(eta) + self.eta * self.i2.eta4(eta))
    }
    fn i2t2d0(&mut self, eta: f64) -> f64 {
        self.eta2 * self.i2.eta1(eta) + self.eta1.powi(2) * self.i2.eta2(eta)
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
            AssocType::Type2B => 2.0 * self.xa.ln() - self.xa + 1.0,
            AssocType::Type3B => {
                2.0 * self.xa.ln() + (2.0 * self.xa - 1.0).ln() - 2.0 * self.xa + 2.0
            }
        }
    }
    fn assoc_t0d1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d1::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t0d1::<1>(self.xa) + self.site_t0d1::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t0d2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d2::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t0d2::<1>(self.xa) + self.site_t0d2::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t0d3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d3::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t0d3::<1>(self.xa) + self.site_t0d3::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t0d4(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t0d4::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t0d4::<1>(self.xa) + self.site_t0d4::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t1d0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d0::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t1d0::<1>(self.xa) + self.site_t1d0::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t1d1(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d1::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t1d1::<1>(self.xa) + self.site_t1d1::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t1d2(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d2::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t1d2::<1>(self.xa) + self.site_t1d2::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t1d3(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t1d3::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t1d3::<1>(self.xa) + self.site_t1d3::<2>(2.0 * self.xa - 1.0)
            }
        }
    }
    fn assoc_t2d0(&mut self) -> f64 {
        match self.assoc_type {
            AssocType::Type0 => 0.0,
            AssocType::Type2B => 2.0 * self.site_t2d0::<1>(self.xa),
            AssocType::Type3B => {
                2.0 * self.site_t2d0::<1>(self.xa) + self.site_t2d0::<2>(2.0 * self.xa - 1.0)
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
    fn xt1(&self) -> f64 {
        (self.rho_num * self.kappa_ab_plus)
            * match self.assoc_type {
                AssocType::Type2B => self.xa.powi(3) / (self.xa - 2.0),
                AssocType::Type3B => {
                    (self.xa * (2.0 * self.xa - 1.0)).powi(2)
                        / (2.0 * self.xa.powi(2) - 4.0 * self.xa + 1.0)
                }
                _ => 0.0,
            }
    }
    fn xt2(&self) -> f64 {
        2.0 * (self.rho_num * self.kappa_ab_plus).powi(2)
            * match self.assoc_type {
                AssocType::Type2B => self.xa.powi(5) / (self.xa - 2.0).powi(3) * (self.xa - 3.0),
                AssocType::Type3B => {
                    (self.xa * (2.0 * self.xa - 1.0)).powi(3)
                        / (2.0 * self.xa.powi(2) - 4.0 * self.xa + 1.0).powi(3)
                        * (4.0 * self.xa.powi(3) - 12.0 * self.xa.powi(2) + 6.0 * self.xa - 1.0)
                }
                _ => 0.0,
            }
    }
    fn xt3(&self) -> f64 {
        6.0 * (self.rho_num * self.kappa_ab_plus).powi(3)
            * match self.assoc_type {
                AssocType::Type2B => {
                    self.xa.powi(7) / (self.xa - 2.0).powi(5)
                        * (self.xa.powi(2) - 6.0 * self.xa + 10.0)
                }
                AssocType::Type3B => {
                    (self.xa * (2.0 * self.xa - 1.0)).powi(4)
                        / (2.0 * self.xa.powi(2) - 4.0 * self.xa + 1.0).powi(5)
                        * (16.0 * self.xa.powi(6) - 96.0 * self.xa.powi(5)
                            + (200.0 * self.xa.powi(4) - 160.0 * self.xa.powi(3))
                            + (62.0 * self.xa.powi(2) - 12.0 * self.xa + 1.0))
                }
                _ => 0.0,
            }
    }
    fn xt4(&self) -> f64 {
        24.0 * (self.rho_num * self.kappa_ab_plus).powi(4)
            * match self.assoc_type {
                AssocType::Type2B => {
                    self.xa.powi(9) / (self.xa - 2.0).powi(7)
                        * (self.xa.powi(3) - 9.0 * self.xa.powi(2) + 29.0 * self.xa - 35.0)
                }
                AssocType::Type3B => {
                    (self.xa * (2.0 * self.xa - 1.0)).powi(5)
                        / (2.0 * self.xa.powi(2) - 4.0 * self.xa + 1.0).powi(7)
                        * (64.0 * self.xa.powi(9) - 576.0 * self.xa.powi(8)
                            + (2080.0 * self.xa.powi(7) - 3808.0 * self.xa.powi(6))
                            + (3696.0 * self.xa.powi(5) - 2084.0 * self.xa.powi(4))
                            + (716.0 * self.xa.powi(3) - 150.0 * self.xa.powi(2))
                            + (18.0 * self.xa - 1.0))
                }
                _ => 0.0,
            }
    }
}
impl PcSaftGlyPure {
    fn x_t0d1(&mut self) -> f64 {
        self.xt1() * self.t_t0d1()
    }
    fn x_t0d2(&mut self) -> f64 {
        self.xt2() * self.t_t0d1().powi(2) + self.xt1() * self.t_t0d2()
    }
    fn x_t0d3(&mut self) -> f64 {
        self.xt3() * self.t_t0d1().powi(3)
            + 3.0 * self.xt2() * self.t_t0d1() * self.t_t0d2()
            + self.xt1() * self.t_t0d3()
    }
    fn x_t0d4(&mut self) -> f64 {
        self.xt4() * self.t_t0d1().powi(4)
            + 6.0 * self.xt3() * self.t_t0d1().powi(2) * self.t_t0d2()
            + 3.0 * self.xt2() * self.t_t0d2().powi(2)
            + 4.0 * self.xt2() * self.t_t0d1() * self.t_t0d3()
            + self.xt1() * self.t_t0d4()
    }
    fn x_t1d0(&mut self) -> f64 {
        self.xt1() * self.t_t1d0()
    }
    fn x_t1d1(&mut self) -> f64 {
        self.xt2() * self.t_t1d0() * self.t_t0d1() + self.xt1() * self.t_t1d1()
    }
    fn x_t1d2(&mut self) -> f64 {
        self.xt3() * self.t_t1d0() * self.t_t0d1().powi(2)
            + 2.0 * self.xt2() * self.t_t1d1() * self.t_t0d1()
            + self.xt2() * self.t_t1d0() * self.t_t0d2()
            + self.xt1() * self.t_t1d2()
    }
    fn x_t1d3(&mut self) -> f64 {
        self.xt4() * self.t_t1d0() * self.t_t0d1().powi(3)
            + 3.0 * self.xt3() * self.t_t1d1() * self.t_t0d1().powi(2)
            + 3.0 * self.xt3() * self.t_t1d0() * self.t_t0d1() * self.t_t0d2()
            + 3.0 * self.xt2() * self.t_t1d2() * self.t_t0d1()
            + 3.0 * self.xt2() * self.t_t1d1() * self.t_t0d2()
            + self.xt2() * self.t_t1d0() * self.t_t0d3()
            + self.xt1() * self.t_t1d3()
    }
    fn x_t2d0(&mut self) -> f64 {
        self.xt2() * self.t_t1d0().powi(2) + self.xt1() * self.t_t2d0()
    }
}
impl PcSaftGlyPure {
    fn t_t0d0(&mut self) -> f64 {
        (self.rho_num * self.kappa_ab_plus)
            * self.gii_gly_t0d0()
            * ((self.epsilon_ab / self.temp).exp() - 1.0)
    }
    fn t_t0d1(&mut self) -> f64 {
        (self.gii_gly_t0d1() + self.gii_gly_t0d0()) * ((self.epsilon_ab / self.temp).exp() - 1.0)
    }
    fn t_t0d2(&mut self) -> f64 {
        (self.gii_gly_t0d2() + 2.0 * self.gii_gly_t0d1())
            * ((self.epsilon_ab / self.temp).exp() - 1.0)
    }
    fn t_t0d3(&mut self) -> f64 {
        (self.gii_gly_t0d3() + 3.0 * self.gii_gly_t0d2())
            * ((self.epsilon_ab / self.temp).exp() - 1.0)
    }
    fn t_t0d4(&self) -> f64 {
        (self.gii_gly_t0d4() + 4.0 * self.gii_gly_t0d3())
            * ((self.epsilon_ab / self.temp).exp() - 1.0)
    }
    fn t_t1d0(&mut self) -> f64 {
        let epsilon_ab_t = self.epsilon_ab / self.temp;
        self.gii_gly_t1d0() * (epsilon_ab_t.exp() - 1.0)
            - self.gii_gly_t0d0() * epsilon_ab_t.exp() * epsilon_ab_t
    }
    fn t_t1d1(&mut self) -> f64 {
        let epsilon_ab_t = self.epsilon_ab / self.temp;
        (self.gii_gly_t1d1() + self.gii_gly_t1d0()) * (epsilon_ab_t.exp() - 1.0)
            - (self.gii_gly_t0d1() + self.gii_gly_t0d0()) * epsilon_ab_t.exp() * epsilon_ab_t
    }
    fn t_t1d2(&mut self) -> f64 {
        let epsilon_ab_t = self.epsilon_ab / self.temp;
        (self.gii_gly_t1d2() + 2.0 * self.gii_gly_t1d1()) * (epsilon_ab_t.exp() - 1.0)
            - (self.gii_gly_t0d2() + 2.0 * self.gii_gly_t0d1()) * epsilon_ab_t.exp() * epsilon_ab_t
    }
    fn t_t1d3(&mut self) -> f64 {
        let epsilon_ab_t = self.epsilon_ab / self.temp;
        (self.gii_gly_t1d3() + 3.0 * self.gii_gly_t1d2()) * (epsilon_ab_t.exp() - 1.0)
            - (self.gii_gly_t0d3() + 3.0 * self.gii_gly_t0d2()) * epsilon_ab_t.exp() * epsilon_ab_t
    }
    fn t_t2d0(&mut self) -> f64 {
        let epsilon_ab_t = self.epsilon_ab / self.temp;
        self.gii_gly_t2d0() * (epsilon_ab_t.exp() - 1.0)
            - 2.0 * self.gii_gly_t1d0() * epsilon_ab_t.exp() * epsilon_ab_t
            + self.gii_gly_t0d0() * epsilon_ab_t.exp() * epsilon_ab_t * (epsilon_ab_t + 2.0)
    }
}
impl PcSaftGlyPure {
    fn gii_gly_t0d0(&mut self) -> f64 {
        (1.0 - self.eta).powi(3).recip()
            * (self.c0 + self.c1 * self.eta + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t0d1(&mut self) -> f64 {
        self.eta / (1.0 - self.eta).powi(4)
            * ((3.0 * self.c0 + self.c1)
                + (2.0 * self.c1 + 2.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t0d2(&mut self) -> f64 {
        self.eta.powi(2) * 2.0 / (1.0 - self.eta).powi(5)
            * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                + (3.0 * self.c1 + 4.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t0d3(&self) -> f64 {
        self.eta.powi(3) * 6.0 / (1.0 - self.eta).powi(6)
            * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                + (4.0 * self.c1 + 6.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t0d4(&self) -> f64 {
        self.eta.powi(4) * 24.0 / (1.0 - self.eta).powi(7)
            * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                + (5.0 * self.c1 + 8.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t1d0(&self) -> f64 {
        self.eta1 / (1.0 - self.eta).powi(4)
            * ((3.0 * self.c0 + self.c1)
                + (2.0 * self.c1 + 2.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t1d1(&self) -> f64 {
        self.eta1 / (1.0 - self.eta).powi(4)
            * ((3.0 * self.c0 + self.c1)
                + (2.0 * self.c1 + 2.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
            + self.eta1 * self.eta * 2.0 / (1.0 - self.eta).powi(5)
                * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                    + (3.0 * self.c1 + 4.0 * self.c2) * self.eta
                    + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t1d2(&self) -> f64 {
        2.0 * self.eta1 * self.eta * 2.0 / (1.0 - self.eta).powi(5)
            * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                + (3.0 * self.c1 + 4.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
            + self.eta1 * self.eta.powi(2) * 6.0 / (1.0 - self.eta).powi(6)
                * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                    + (4.0 * self.c1 + 6.0 * self.c2) * self.eta
                    + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t1d3(&self) -> f64 {
        3.0 * self.eta1 * self.eta.powi(2) * 6.0 / (1.0 - self.eta).powi(6)
            * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                + (4.0 * self.c1 + 6.0 * self.c2) * self.eta
                + self.c2 * self.eta.powi(2))
            + self.eta1 * self.eta.powi(3) * 24.0 / (1.0 - self.eta).powi(7)
                * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                    + (5.0 * self.c1 + 8.0 * self.c2) * self.eta
                    + self.c2 * self.eta.powi(2))
    }
    fn gii_gly_t2d0(&self) -> f64 {
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
impl PcSaftGlyPure {
    pub fn check_derivatives(&mut self, print_val: bool) {
        let (t, d) = (self.temp, self.rho_num);
        if print_val {
            println!("[rT0D0 == rT0D0] calc_rT0D0() ={}", self.r_t0d0(t, d));
        }
        let compare_val = |val_calc: f64, val_diff: f64| {
            assert_eq!(
                &val_calc.abs().to_string()[0..11],
                &val_diff.abs().to_string()[0..11]
            )
        };
        // derivative for density
        let val_calc = self.r_t0d1(t, d) / d;
        let val_diff = romberg_diff(|dx: f64| self.r_t0d0(t, dx), d);
        if print_val {
            println!("[rT0D1 == rT0D1] calc_rT0D1() ={}", val_calc);
            println!("[rT0D0 -> rT0D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density
        let val_calc = self.r_t0d2(t, d) / d.powi(2);
        let val_diff = romberg_diff(|dx: f64| self.r_t0d1(t, dx) / dx, d);
        if print_val {
            println!("[rT0D2 == rT0D2] calc_rT0D2() ={}", val_calc);
            println!("[rT0D1 -> rT0D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density
        let val_calc = self.r_t0d3(t, d) / d.powi(3);
        let val_diff = romberg_diff(|dx: f64| self.r_t0d2(t, dx) / dx.powi(2), d);
        if print_val {
            println!("[rT0D3 == rT0D3] calc_rT0D3() ={}", val_calc);
            println!("[rT0D2 -> rT0D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density+density+density
        let val_calc = self.r_t0d4(t, d) / d.powi(4);
        let val_diff = romberg_diff(|dx: f64| self.r_t0d3(t, dx) / dx.powi(3), d);
        if print_val {
            println!("[rT0D4 == rT0D4] calc_rT0D4() ={}", val_calc);
            println!("[rT0D3 -> rT0D4] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature
        let val_calc = self.r_t1d0(t, d) / t;
        let val_diff = romberg_diff(|tx: f64| self.r_t0d0(tx, d), t);
        if print_val {
            println!("[rT1D0 == rT1D0] calc_rT1D0() ={}", val_calc);
            println!("[rT0D0 -> rT1D0] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density
        let val_calc = self.r_t1d1(t, d) / t / d;
        let val_diff = romberg_diff(|dx: f64| self.r_t1d0(t, dx) / t, d);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_rT1D1() ={}", val_calc);
            println!("[rT1D0 -> rT1D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.r_t0d1(tx, d) / d, t);
        if print_val {
            println!("[rT1D1 == rT1D1] calc_rT1D1() ={}", val_calc);
            println!("[rT0D1 -> rT1D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density
        let val_calc = self.r_t1d2(t, d) / t / d.powi(2);
        let val_diff = romberg_diff(|dx: f64| self.r_t1d1(t, dx) / t / dx, d);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_rT1D2() ={}", val_calc);
            println!("[rT1D1 -> rT1D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.r_t0d2(tx, d) / d.powi(2), t);
        if print_val {
            println!("[rT1D2 == rT1D2] calc_rT1D2() ={}", val_calc);
            println!("[rT0D2 -> rT1D2] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+density+density+density
        let val_calc = self.r_t1d3(t, d) / t / d.powi(3);
        let val_diff = romberg_diff(|dx: f64| self.r_t1d2(t, dx) / t / dx.powi(2), d);
        if print_val {
            println!("[rT1D3 == rT1D3] calc_rT1D3() ={}", val_calc);
            println!("[rT1D2 -> rT1D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        let val_diff = romberg_diff(|tx: f64| self.r_t0d3(tx, d) / d.powi(3), t);
        if print_val {
            println!("[rT1D3 == rT1D3] calc_rT1D3() ={}", val_calc);
            println!("[rT0D3 -> rT1D3] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for temperature+temperature
        let val_calc = self.r_t2d0(t, d) / t.powi(2);
        let val_diff = romberg_diff(|tx: f64| self.r_t1d0(tx, d) / tx, t);
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
        let temp_min = (0.6 * fluid.t().unwrap()).floor() as u32;
        let temp_max = fluid.t().unwrap().ceil() as u32;
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
        methanol.set_2b_assoc_type(0.035176, 2899.5);
        methanol.set_gly_params(1.0, 1.0, 1.0);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
        let mut methanol = PcSaftGlyPure::new_fluid(1.5255, 3.23, 188.9);
        methanol.set_3b_assoc_type(0.035176, 2899.5);
        methanol.set_gly_params(1.0, 1.0, 1.0);
        methanol.td_unchecked(300.0, 24514.0);
        methanol.check_derivatives(false);
    }
}
