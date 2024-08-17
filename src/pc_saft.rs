use thiserror::Error;
#[derive(Debug, Error)]
enum PcSaftPureErr {
    #[error("c_flash diverge")]
    NotConvForC,
    #[error("tp_flash diverge")]
    NotConvForTP,
    #[error("property not in single phase")]
    NotInSinglePhase,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
}
use anyhow::anyhow;
use pyo3::{pyclass, pymethods};
/// PC-SAFT EOS
/// ```
/// ```
#[pyclass]
#[allow(non_snake_case)]
pub struct PcSaftPure {
    m: f64,
    sigma: f64,
    epsilon: f64,
    T: f64,
    rho_num: f64,
    rhov_num: f64,
    rhol_num: f64,
    is_single_phase: bool,
    m1: f64,  // (m-1)/m
    m12: f64, // (m-1)/m * (m-2)/m
    // changed from T and rho_num
    eta: f64,
    eta1: f64,
    m2e1s3: f64,
    m2e2s3: f64,
}
#[pymethods]
#[allow(non_snake_case)]
impl PcSaftPure {
    #[new]
    pub fn new_fluid(m: f64, sigma: f64, epsilon: f64) -> Self {
        Self {
            m,
            sigma,
            epsilon,
            T: 1.0,
            rho_num: 1E-10,
            rhov_num: 0.0,
            rhol_num: 0.0,
            is_single_phase: true,
            m1: (m - 1.0) / m,                      // (m-1)/m
            m12: (m - 1.0) * (m - 2.0) / m.powi(2), // (m-1)/m * (m-2)/m
            // changed from T and rho_num
            eta: 0.0,
            eta1: 0.0,
            m2e1s3: 0.0,
            m2e2s3: 0.0,
        }
    }
    pub fn check_derivatives(&mut self) {
        let T = self.T;
        let rho = self.rho_num;
        // derivative for density
        let rT0D0_rT0D1 = |rhox: f64| self.calc_rT0D0(T, rhox);
        println!(
            "[rT0D0 -> rT0D1] romberg_diff ={}",
            romberg_diff(rT0D0_rT0D1, rho)
        );
        println!(
            "[rT0D1 == rT0D1] calc_rT0D1() ={}",
            self.calc_rT0D1(T, rho) / rho
        );
        // derivative for density+density
        let rT0D1_rT0D2 = |rhox: f64| self.calc_rT0D1(T, rhox) / rhox;
        println!(
            "[rT0D1 -> rT0D2] romberg_diff ={}",
            romberg_diff(rT0D1_rT0D2, rho)
        );
        println!(
            "[rT0D2 == rT0D2] calc_rT0D2() ={}",
            self.calc_rT0D2(T, rho) / rho.powi(2)
        );
        // derivative for density+density+density
        let rT0D2_rT0D3 = |rhox: f64| self.calc_rT0D2(T, rhox) / rhox.powi(2);
        println!(
            "[rT0D2 -> rT0D3] romberg_diff ={}",
            romberg_diff(rT0D2_rT0D3, rho)
        );
        println!(
            "[rT0D3 == rT0D3] calc_rT0D3() ={}",
            self.calc_rT0D3(T, rho) / rho.powi(3)
        );
        // derivative for density+density+density+density
        let rT0D3_rT0D4 = |rhox: f64| self.calc_rT0D3(T, rhox) / rhox.powi(3);
        println!(
            "[rT0D3 -> rT0D4] romberg_diff ={}",
            romberg_diff(rT0D3_rT0D4, rho)
        );
        println!(
            "[rT0D4 == rT0D4] calc_rT0D4() ={}",
            self.calc_rT0D4(T, rho) / rho.powi(4)
        );
        // derivative for temperature+density
        let rT0D1_rT1D1 = |Tx: f64| self.calc_rT0D1(Tx, rho) / rho;
        println!(
            "[rT0D1 -> rT1D1] romberg_diff ={}",
            romberg_diff(rT0D1_rT1D1, T)
        );
        println!(
            "[rT1D1 == rT1D1] calc_rT1D1() ={}",
            self.calc_rT1D1(T, rho) / T / rho
        );
        // derivative for temperature+density+density
        let rT0D2_rT1D2 = |Tx: f64| self.calc_rT0D2(Tx, rho) / rho.powi(2);
        println!(
            "[rT0D2 -> rT1D2] romberg_diff ={}",
            romberg_diff(rT0D2_rT1D2, T)
        );
        println!(
            "[rT1D2 == rT1D2] calc_rT1D2() ={}",
            self.calc_rT1D2(T, rho) / T / rho.powi(2)
        );
        let rT1D1_rT1D2 = |rhox: f64| self.calc_rT1D1(T, rhox) / T / rhox;
        println!(
            "[rT1D1 -> rT1D2] romberg_diff ={}",
            romberg_diff(rT1D1_rT1D2, rho)
        );
        println!(
            "[rT1D2 == rT1D2] calc_rT1D2() ={}",
            self.calc_rT1D2(T, rho) / T / rho.powi(2)
        );
        // derivative for temperature+density+density+density
        let rT0D3_rT1D3 = |Tx: f64| self.calc_rT0D3(Tx, rho) / rho.powi(3);
        println!(
            "[rT0D3 -> rT1D3] romberg_diff ={}",
            romberg_diff(rT0D3_rT1D3, T)
        );
        println!(
            "[rT1D3 == rT1D3] calc_rT1D3() ={}",
            self.calc_rT1D3(T, rho) / T / rho.powi(3)
        );
        let rT1D2_rT1D3 = |rhox: f64| self.calc_rT1D2(T, rhox) / T / rhox.powi(2);
        println!(
            "[rT1D2 -> rT1D3] romberg_diff ={}",
            romberg_diff(rT1D2_rT1D3, rho)
        );
        println!(
            "[rT1D3 == rT1D3] calc_rT1D3() ={}",
            self.calc_rT1D3(T, rho) / T / rho.powi(3)
        );
    }
    pub fn td_unchecked(&mut self, T: f64, rho_mol: f64) {
        self.set_temperature_and_number_density(T, rho_mol * NA / 1E30);
        self.is_single_phase = true;
    }
    pub fn c_flash(&mut self) -> anyhow::Result<()> {
        // Iteration from T_c = 1000 eta_c = 1E-10
        let mut Tc = 1000.0;
        let mut rhoc = (6E-10 / PI / self.m)
            / (self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / Tc).exp())).powi(3);
        // Define variables
        let mut Dp_Drho_T = self.calc_Dp_Drho_T(Tc, rhoc);
        let mut D2p_DTrho;
        let mut D2p_Drho2_T = self.calc_D2p_Drho2_T(Tc, rhoc);
        let mut D3p_DTrho2;
        let mut D3p_Drho3_T;
        for _i in 1..100000 {
            D2p_DTrho = self.calc_D2p_DTrho(Tc, rhoc);
            D3p_DTrho2 = self.calc_D3p_DTrho2(Tc, rhoc);
            D3p_Drho3_T = self.calc_D3p_Drho3_T(Tc, rhoc);
            Tc -= (Dp_Drho_T * D3p_Drho3_T - D2p_Drho2_T * D2p_Drho2_T)
                / (D2p_DTrho * D3p_Drho3_T - D3p_DTrho2 * D2p_Drho2_T);
            rhoc -= (Dp_Drho_T * D3p_DTrho2 - D2p_Drho2_T * D2p_DTrho)
                / (D2p_Drho2_T * D3p_DTrho2 - D3p_Drho3_T * D2p_DTrho);
            Dp_Drho_T = self.calc_Dp_Drho_T(Tc, rhoc);
            D2p_Drho2_T = self.calc_D2p_Drho2_T(Tc, rhoc);
            if Dp_Drho_T.abs() < 1E3 && D2p_Drho2_T.abs() < 1E6 {
                self.set_temperature_and_number_density(Tc, rhoc);
                self.is_single_phase = true;
                return Ok(());
            }
        }
        Err(anyhow!(PcSaftPureErr::NotConvForC))
    }
    pub fn tp_flash(&mut self, T: f64, p: f64) -> anyhow::Result<()> {
        let d3 = (self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / T).exp())).powi(3);
        let mut p_diff: f64;
        // Iteration from gas phase: eta = 1E-10
        let mut rhov_num = 6.0 * 1E-10 / PI / self.m / d3;
        let mut Dp_Drhov_T = -1.0;
        loop {
            p_diff = self.calc_p(T, rhov_num) - p;
            if (p_diff / p).abs() < 1E-9 {
                break;
            }
            Dp_Drhov_T = self.calc_Dp_Drho_T(T, rhov_num);
            if Dp_Drhov_T.is_sign_negative() {
                break;
            }
            rhov_num -= p_diff / Dp_Drhov_T;
            if rhov_num.is_sign_negative() {
                break;
            }
        }
        let lnphi_v = if rhov_num.is_sign_negative() || Dp_Drhov_T.is_sign_negative() {
            1E16
        } else {
            self.calc_lnphi(T, rhov_num)
        };
        // Iteration from liquid phase: eta = 0.5
        let mut rhol_num = 6.0 * 0.5 / PI / self.m / d3;
        let mut Dp_Drhol_T = -1.0;
        loop {
            p_diff = self.calc_p(T, rhol_num) - p;
            if (p_diff / p).abs() < 1E-9 {
                break;
            }
            Dp_Drhol_T = self.calc_Dp_Drho_T(T, rhol_num);
            if Dp_Drhol_T.is_sign_negative() {
                break;
            }
            rhol_num -= p_diff / Dp_Drhol_T;
            if rhol_num.is_sign_negative() {
                break;
            }
        }
        let lnphi_l = if rhol_num.is_sign_negative() || Dp_Drhol_T.is_sign_negative() {
            1E16
        } else {
            self.calc_lnphi(T, rhol_num)
        };
        // Select the correct output
        if lnphi_v < 1E16 && lnphi_l < 1E16 {
            if lnphi_v < lnphi_l {
                self.set_temperature_and_number_density(T, rhov_num);
            } else {
                self.set_temperature_and_number_density(T, rhol_num);
            }
            self.is_single_phase = true;
            Ok(())
        } else if lnphi_v < 1E16 && lnphi_l == 1E16 {
            self.set_temperature_and_number_density(T, rhov_num);
            self.is_single_phase = true;
            Ok(())
        } else if lnphi_l < 1E16 && lnphi_v == 1E16 {
            self.set_temperature_and_number_density(T, rhol_num);
            self.is_single_phase = true;
            Ok(())
        } else {
            Err(anyhow!(PcSaftPureErr::NotConvForTP))
        }
    }
    pub fn T(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.T)
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn p(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.calc_p(self.T, self.rho_num))
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn rho(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.rho_num * 1E30 / NA)
        } else {
            Err(anyhow!(PcSaftPureErr::OnlyInSinglePhase))
        }
    }
    pub fn T_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::NotInSinglePhase))
        } else {
            Ok(self.T)
        }
    }
    pub fn p_s(&mut self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::NotInSinglePhase))
        } else {
            Ok((self.calc_p(self.T, self.rhov_num) + self.calc_p(self.T, self.rhol_num)) / 2.0)
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::NotInSinglePhase))
        } else {
            Ok(self.rhov_num * 1E30 / NA)
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PcSaftPureErr::NotInSinglePhase))
        } else {
            Ok(self.rhol_num * 1E30 / NA)
        }
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn calc_p(&mut self, T: f64, rho_num: f64) -> f64 {
        (1E30 / NA * R * T) * self.rho_num * (1.0 + self.calc_rT0D1(T, rho_num))
    }
    fn calc_Dp_Drho_T(&mut self, T: f64, rho_num: f64) -> f64 {
        (1E30 / NA * R * T)
            * (1.0 + 2.0 * self.calc_rT0D1(T, rho_num) + self.calc_rT0D2(T, rho_num))
    }
    fn calc_D2p_Drho2_T(&mut self, T: f64, rho_num: f64) -> f64 {
        (1E30 / NA * R * T)
            * (2.0 * self.calc_rT0D1(T, rho_num)
                + 4.0 * self.calc_rT0D2(T, rho_num)
                + self.calc_rT0D3(T, rho_num))
            / rho_num
    }
    fn calc_D3p_Drho3_T(&mut self, T: f64, rho_num: f64) -> f64 {
        (1E30 / NA * R * T)
            * (6.0 * self.calc_rT0D2(T, rho_num)
                + 6.0 * self.calc_rT0D3(T, rho_num)
                + self.calc_rT0D4(T, rho_num))
            / rho_num.powi(2)
    }
    fn calc_D2p_DTrho(&mut self, T: f64, rho_num: f64) -> f64 {
        (1E30 / NA * R * T)
            * (1.0
                + 2.0 * self.calc_rT0D1(T, rho_num)
                + self.calc_rT0D2(T, rho_num)
                + 2.0 * self.calc_rT1D1(T, rho_num)
                + self.calc_rT1D2(T, rho_num))
    }
    fn calc_D3p_DTrho2(&mut self, T: f64, rho_num: f64) -> f64 {
        (1E30 / NA * R * T)
            * (2.0 * self.calc_rT0D1(T, rho_num)
                + 4.0 * self.calc_rT0D2(T, rho_num)
                + self.calc_rT0D3(T, rho_num)
                + 2.0 * self.calc_rT1D1(T, rho_num)
                + 4.0 * self.calc_rT1D2(T, rho_num)
                + self.calc_rT1D3(T, rho_num))
            / rho_num
    }
    fn calc_lnphi(&mut self, T: f64, rho_num: f64) -> f64 {
        self.calc_rT0D0(T, rho_num) + self.calc_rT0D1(T, rho_num)
            - (1.0 + self.calc_rT0D1(T, rho_num)).ln()
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn set_temperature_and_number_density(&mut self, T: f64, rho_num: f64) {
        if T != self.T {
            self.T = T;
            self.m2e1s3 = self.m.powi(2) * (self.epsilon / T) * self.sigma.powi(3);
            self.m2e2s3 = self.m.powi(2) * (self.epsilon / T).powi(2) * self.sigma.powi(3);
            self.rho_num = rho_num;
            let d = self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / T).exp());
            let d1 = -0.36 * self.sigma * (-3.0 * self.epsilon / T).exp() * self.epsilon / T;
            self.eta = PI / 6.0 * rho_num * self.m * d.powi(3);
            self.eta1 = PI / 2.0 * rho_num * self.m * d.powi(2) * d1;
        } else if rho_num != self.rho_num {
            self.rho_num = rho_num;
            let d = self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / T).exp());
            let d1 = -0.36 * self.sigma * (-3.0 * self.epsilon / T).exp() * self.epsilon / T;
            self.eta = PI / 6.0 * rho_num * self.m * d.powi(3);
            self.eta1 = PI / 2.0 * rho_num * self.m * d.powi(2) * d1;
        }
    }
    fn calc_rT0D0(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT0D0(self.eta) + self.dispT0D0(self.eta)
    }
    fn calc_rT0D1(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT0D1(self.eta) + self.dispT0D1(self.eta)
    }
    fn calc_rT0D2(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT0D2(self.eta) + self.dispT0D2(self.eta)
    }
    fn calc_rT0D3(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT0D3(self.eta) + self.dispT0D3(self.eta)
    }
    fn calc_rT0D4(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT0D4(self.eta) + self.dispT0D4(self.eta)
    }
    fn calc_rT1D1(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT1D1(self.eta) + self.dispT1D1(self.eta)
    }
    fn calc_rT1D2(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT1D2(self.eta) + self.dispT1D2(self.eta)
    }
    fn calc_rT1D3(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.hcT1D3(self.eta) + self.dispT1D3(self.eta)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn hcT0D0(&self, eta: f64) -> f64 {
        self.m * self.hsT0D0(eta) + (1.0 - self.m) * self.giiT0D0(eta).ln()
    }
    fn hcT0D1(&self, eta: f64) -> f64 {
        eta * (self.m * self.hsT0D1(eta) + (1.0 - self.m) / self.giiT0D0(eta) * self.giiT0D1(eta))
    }
    fn hcT0D2(&self, eta: f64) -> f64 {
        let giiT0D0: f64 = self.giiT0D0(eta);
        eta.powi(2)
            * (self.m * self.hsT0D2(eta)
                + (1.0 - self.m)
                    * (-1.0 / giiT0D0.powi(2) * self.giiT0D1(eta).powi(2)
                        + 1.0 / giiT0D0 * self.giiT0D2(eta)))
    }
    fn hcT0D3(&self, eta: f64) -> f64 {
        let giiT0D0: f64 = self.giiT0D0(eta);
        let giiT0D1: f64 = self.giiT0D1(eta);
        eta.powi(3)
            * (self.m * self.hsT0D3(eta)
                + (1.0 - self.m)
                    * (2.0 / giiT0D0.powi(3) * giiT0D1.powi(3)
                        - 3.0 / giiT0D0.powi(2) * giiT0D1 * self.giiT0D2(eta)
                        + 1.0 / giiT0D0 * self.giiT0D3(eta)))
    }
    fn hcT0D4(&self, eta: f64) -> f64 {
        let giiT0D0: f64 = self.giiT0D0(eta);
        let giiT0D1: f64 = self.giiT0D1(eta);
        let giiT0D2: f64 = self.giiT0D2(eta);
        eta.powi(4)
            * (self.m * self.hsT0D4(eta)
                + (1.0 - self.m)
                    * (-6.0 / giiT0D0.powi(4) * giiT0D1.powi(4)
                        + 12.0 / giiT0D0.powi(3) * giiT0D1.powi(2) * giiT0D2
                        - 4.0 / giiT0D0.powi(2) * giiT0D1 * self.giiT0D3(eta)
                        - 3.0 / giiT0D0.powi(2) * giiT0D2.powi(2)
                        + 1.0 / giiT0D0 * self.giiT0D4(eta)))
    }
    fn hcT1D1(&self, eta: f64) -> f64 {
        let giiT0D0: f64 = self.giiT0D0(eta);
        self.eta1
            * eta
            * (self.m * self.hsT1D1(eta)
                + (1.0 - self.m)
                    * (-1.0 / giiT0D0.powi(2) * self.giiT1D0(eta) * self.giiT0D1(eta)
                        + 1.0 / giiT0D0 * self.giiT1D1(eta)))
    }
    fn hcT1D2(&self, eta: f64) -> f64 {
        let giiT0D0: f64 = self.giiT0D0(eta);
        let giiT0D1: f64 = self.giiT0D1(eta);
        let giiT1D0: f64 = self.giiT1D0(eta);
        self.eta1
            * eta.powi(2)
            * (self.m * self.hsT1D2(eta)
                + (1.0 - self.m)
                    * (2.0 / giiT0D0.powi(3) * giiT1D0 * giiT0D1.powi(2)
                        - 2.0 / giiT0D0.powi(2) * self.giiT1D1(eta) * giiT0D1
                        - 1.0 / giiT0D0.powi(2) * giiT1D0 * self.giiT0D2(eta)
                        + 1.0 / giiT0D0 * self.giiT1D2(eta)))
    }
    fn hcT1D3(&self, eta: f64) -> f64 {
        let giiT0D0: f64 = self.giiT0D0(eta);
        let giiT0D1: f64 = self.giiT0D1(eta);
        let giiT0D2: f64 = self.giiT0D2(eta);
        let giiT1D0: f64 = self.giiT1D0(eta);
        let giiT1D1: f64 = self.giiT1D1(eta);
        self.eta1
            * eta.powi(3)
            * (self.m * self.hsT1D3(eta)
                + (1.0 - self.m)
                    * (-6.0 / giiT0D0.powi(4) * giiT1D0 * giiT0D1.powi(3)
                        + 6.0 / giiT0D0.powi(3) * giiT1D1 * giiT0D1.powi(2)
                        + 6.0 / giiT0D0.powi(3) * giiT1D0 * giiT0D1 * giiT0D2
                        - 3.0 / giiT0D0.powi(2) * giiT1D1 * giiT0D2
                        - 3.0 / giiT0D0.powi(2) * self.giiT1D2(eta) * giiT0D1
                        - 1.0 / giiT0D0.powi(2) * giiT1D0 * self.giiT0D3(eta)
                        + 1.0 / giiT0D0 * self.giiT1D3(eta)))
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn hsT0D0(&self, eta: f64) -> f64 {
        eta * (4.0 - 3.0 * eta) / (1.0 - eta).powi(2)
    }
    fn hsT0D1(&self, eta: f64) -> f64 {
        (4.0 - 2.0 * eta) / (1.0 - eta).powi(3)
    }
    fn hsT0D2(&self, eta: f64) -> f64 {
        (10.0 - 4.0 * eta) / (1.0 - eta).powi(4)
    }
    fn hsT0D3(&self, eta: f64) -> f64 {
        (36.0 - 12.0 * eta) / (1.0 - eta).powi(5)
    }
    fn hsT0D4(&self, eta: f64) -> f64 {
        (168.0 - 48.0 * eta) / (1.0 - eta).powi(6)
    }
    fn hsT1D1(&self, eta: f64) -> f64 {
        (4.0 - 2.0 * eta) / (1.0 - eta).powi(3) / eta + (10.0 - 4.0 * eta) / (1.0 - eta).powi(4)
    }
    fn hsT1D2(&self, eta: f64) -> f64 {
        2.0 * (10.0 - 4.0 * eta) / (1.0 - eta).powi(4) / eta
            + (36.0 - 12.0 * eta) / (1.0 - eta).powi(5)
    }
    fn hsT1D3(&self, eta: f64) -> f64 {
        3.0 * (36.0 - 12.0 * eta) / (1.0 - eta).powi(5) / eta
            + (168.0 - 48.0 * eta) / (1.0 - eta).powi(6)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn giiT0D0(&self, eta: f64) -> f64 {
        (1.0 - 0.5 * eta) / (1.0 - eta).powi(3)
    }
    fn giiT0D1(&self, eta: f64) -> f64 {
        (2.5 - eta) / (1.0 - eta).powi(4)
    }
    fn giiT0D2(&self, eta: f64) -> f64 {
        (9.0 - 3.0 * eta) / (1.0 - eta).powi(5)
    }
    fn giiT0D3(&self, eta: f64) -> f64 {
        (42.0 - 12.0 * eta) / (1.0 - eta).powi(6)
    }
    fn giiT0D4(&self, eta: f64) -> f64 {
        (240.0 - 60.0 * eta) / (1.0 - eta).powi(7)
    }
    fn giiT1D0(&self, eta: f64) -> f64 {
        (2.5 - eta) / (1.0 - eta).powi(4)
    }
    fn giiT1D1(&self, eta: f64) -> f64 {
        (2.5 - eta) / (1.0 - eta).powi(4) / eta + (9.0 - 3.0 * eta) / (1.0 - eta).powi(5)
    }
    fn giiT1D2(&self, eta: f64) -> f64 {
        2.0 * (9.0 - 3.0 * eta) / (1.0 - eta).powi(5) / eta
            + (42.0 - 12.0 * eta) / (1.0 - eta).powi(6)
    }
    fn giiT1D3(&self, eta: f64) -> f64 {
        3.0 * (42.0 - 12.0 * eta) / (1.0 - eta).powi(6) / eta
            + (240.0 - 60.0 * eta) / (1.0 - eta).powi(7)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn dispT0D0(&self, eta: f64) -> f64 {
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * self.i1T0D0(eta)
                + self.m * self.m2e2s3 * self.c1T0D0(eta) * self.i2T0D0(eta))
    }
    fn dispT0D1(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let i2T0D0 = self.i2T0D0(eta);
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (self.i1T0D0(eta) + self.i1T0D1(eta))
                + self.m
                    * self.m2e2s3
                    * (c1T0D0 * i2T0D0 + self.c1T0D1(eta) * i2T0D0 + c1T0D0 * self.i2T0D1(eta)))
    }
    fn dispT0D2(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let c1T0D1 = self.c1T0D1(eta);
        let i2T0D0 = self.i2T0D0(eta);
        let i2T0D1 = self.i2T0D1(eta);
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (2.0 * self.i1T0D1(eta) + self.i1T0D2(eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * c1T0D1 * i2T0D0
                        + 2.0 * c1T0D0 * i2T0D1
                        + self.c1T0D2(eta) * i2T0D0
                        + 2.0 * c1T0D1 * i2T0D1
                        + c1T0D0 * self.i2T0D2(eta)))
    }
    fn dispT0D3(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let c1T0D1 = self.c1T0D1(eta);
        let c1T0D2 = self.c1T0D2(eta);
        let i2T0D0 = self.i2T0D0(eta);
        let i2T0D1 = self.i2T0D1(eta);
        let i2T0D2 = self.i2T0D2(eta);
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (3.0 * self.i1T0D2(eta) + self.i1T0D3(eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * c1T0D2 * i2T0D0
                        + 6.0 * c1T0D1 * i2T0D1
                        + 3.0 * c1T0D0 * i2T0D2
                        + self.c1T0D3(eta) * i2T0D0
                        + 3.0 * c1T0D2 * i2T0D1
                        + 3.0 * c1T0D1 * i2T0D2
                        + c1T0D0 * self.i2T0D3(eta)))
    }
    fn dispT0D4(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let c1T0D1 = self.c1T0D1(eta);
        let c1T0D2 = self.c1T0D2(eta);
        let c1T0D3 = self.c1T0D3(eta);
        let i2T0D0 = self.i2T0D0(eta);
        let i2T0D1 = self.i2T0D1(eta);
        let i2T0D2 = self.i2T0D2(eta);
        let i2T0D3 = self.i2T0D3(eta);
        -PI * self.rho_num
            * (2.0 * self.m2e1s3 * (4.0 * self.i1T0D3(eta) + self.i1T0D4(eta))
                + self.m
                    * self.m2e2s3
                    * (4.0 * c1T0D3 * i2T0D0
                        + 12.0 * c1T0D2 * i2T0D1
                        + 12.0 * c1T0D1 * i2T0D2
                        + 4.0 * c1T0D0 * i2T0D3
                        + self.c1T0D4(eta) * i2T0D0
                        + 4.0 * c1T0D3 * i2T0D1
                        + 6.0 * c1T0D2 * i2T0D2
                        + 4.0 * c1T0D1 * i2T0D3
                        + c1T0D0 * self.i2T0D4(eta)))
    }
    fn dispT1D1(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let c1T0D1 = self.c1T0D1(eta);
        let c1T1D0 = self.c1T1D0(eta);
        let i2T0D0 = self.i2T0D0(eta);
        let i2T0D1 = self.i2T0D1(eta);
        let i2T1D0 = self.i2T1D0(eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (self.i1T1D0(eta) - self.i1T0D0(eta) + self.i1T1D1(eta) - self.i1T0D1(eta))
                + self.m
                    * self.m2e2s3
                    * (c1T1D0 * i2T0D0 + c1T0D0 * i2T1D0 - 2.0 * c1T0D0 * i2T0D0
                        + self.c1T1D1(eta) * i2T0D0
                        + c1T0D1 * i2T1D0
                        - 2.0 * c1T0D1 * i2T0D0
                        + c1T1D0 * i2T0D1
                        + c1T0D0 * self.i2T1D1(eta)
                        - 2.0 * c1T0D0 * i2T0D1))
    }
    fn dispT1D2(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let c1T0D1 = self.c1T0D1(eta);
        let c1T0D2 = self.c1T0D2(eta);
        let c1T1D0 = self.c1T1D0(eta);
        let c1T1D1 = self.c1T1D1(eta);
        let i2T0D0 = self.i2T0D0(eta);
        let i2T0D1 = self.i2T0D1(eta);
        let i2T0D2 = self.i2T0D2(eta);
        let i2T1D0 = self.i2T1D0(eta);
        let i2T1D1 = self.i2T1D1(eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (2.0 * self.i1T1D1(eta) - 2.0 * self.i1T0D1(eta) + self.i1T1D2(eta)
                    - self.i1T0D2(eta))
                + self.m
                    * self.m2e2s3
                    * (2.0 * c1T1D1 * i2T0D0 + 2.0 * c1T0D1 * i2T1D0 - 4.0 * c1T0D1 * i2T0D0
                        + 2.0 * c1T1D0 * i2T0D1
                        + 2.0 * c1T0D0 * i2T1D1
                        - 4.0 * c1T0D0 * i2T0D1
                        + self.c1T1D2(eta) * i2T0D0
                        + c1T0D2 * i2T1D0
                        - 2.0 * c1T0D2 * i2T0D0
                        + 2.0 * c1T1D1 * i2T0D1
                        + 2.0 * c1T0D1 * i2T1D1
                        - 4.0 * c1T0D1 * i2T0D1
                        + c1T1D0 * i2T0D2
                        + c1T0D0 * self.i2T1D2(eta)
                        - 2.0 * c1T0D0 * i2T0D2))
    }
    fn dispT1D3(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let c1T0D1 = self.c1T0D1(eta);
        let c1T0D2 = self.c1T0D2(eta);
        let c1T0D3 = self.c1T0D3(eta);
        let c1T1D0 = self.c1T1D0(eta);
        let c1T1D1 = self.c1T1D1(eta);
        let c1T1D2 = self.c1T1D2(eta);
        let i2T0D0 = self.i2T0D0(eta);
        let i2T0D1 = self.i2T0D1(eta);
        let i2T0D2 = self.i2T0D2(eta);
        let i2T0D3 = self.i2T0D3(eta);
        let i2T1D0 = self.i2T1D0(eta);
        let i2T1D1 = self.i2T1D1(eta);
        let i2T1D2 = self.i2T1D2(eta);
        -PI * self.rho_num
            * (2.0
                * self.m2e1s3
                * (3.0 * self.i1T1D2(eta) - 3.0 * self.i1T0D2(eta) + self.i1T1D3(eta)
                    - self.i1T0D3(eta))
                + self.m
                    * self.m2e2s3
                    * (3.0 * self.c1T1D2(eta) * i2T0D0 + 3.0 * c1T0D2 * i2T1D0
                        - 6.0 * c1T0D2 * i2T0D0
                        + 6.0 * c1T1D1 * i2T0D1
                        + 6.0 * c1T0D1 * i2T1D1
                        - 12.0 * c1T0D1 * i2T0D1
                        + 3.0 * c1T1D0 * i2T0D2
                        + 3.0 * c1T0D0 * self.i2T1D2(eta)
                        - 6.0 * c1T0D0 * i2T0D2
                        + self.c1T1D3(eta) * i2T0D0
                        + c1T0D3 * i2T1D0
                        - 2.0 * c1T0D3 * i2T0D0
                        + 3.0 * c1T1D2 * i2T0D1
                        + 3.0 * c1T0D2 * i2T1D1
                        - 6.0 * c1T0D2 * i2T0D1
                        + 3.0 * c1T1D1 * i2T0D2
                        + 3.0 * c1T0D1 * i2T1D2
                        - 6.0 * c1T0D1 * i2T0D2
                        + c1T1D0 * i2T0D3
                        + c1T0D0 * self.i2T1D3(eta)
                        - 2.0 * c1T0D0 * i2T0D3))
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn c1T0D0(&self, eta: f64) -> f64 {
        1.0 / self.cT0D0(eta)
    }
    fn c1T0D1(&self, eta: f64) -> f64 {
        -self.cT0D1(eta) / self.cT0D0(eta).powi(2)
    }
    fn c1T0D2(&self, eta: f64) -> f64 {
        let cT0D0 = self.cT0D0(eta);
        2.0 * self.cT0D1(eta).powi(2) / cT0D0.powi(3) - self.cT0D2(eta) / cT0D0.powi(2)
    }
    fn c1T0D3(&self, eta: f64) -> f64 {
        let cT0D0 = self.cT0D0(eta);
        let cT0D1 = self.cT0D1(eta);
        -6.0 * cT0D1.powi(3) / cT0D0.powi(4) + 6.0 * cT0D1 * self.cT0D2(eta) / cT0D0.powi(3)
            - self.cT0D3(eta) / cT0D0.powi(2)
    }
    fn c1T0D4(&self, eta: f64) -> f64 {
        let cT0D0 = self.cT0D0(eta);
        let cT0D1 = self.cT0D1(eta);
        let cT0D2 = self.cT0D2(eta);
        24.0 * cT0D1.powi(4) / cT0D0.powi(5) - 36.0 * cT0D1.powi(2) * cT0D2 / cT0D0.powi(4)
            + 8.0 * cT0D1 * self.cT0D3(eta) / cT0D0.powi(3)
            + 6.0 * cT0D2.powi(2) / cT0D0.powi(3)
            - self.cT0D4(eta) / cT0D0.powi(2)
    }
    fn c1T1D0(&self, eta: f64) -> f64 {
        -self.cT1D0(eta) / self.cT0D0(eta).powi(2)
    }
    fn c1T1D1(&self, eta: f64) -> f64 {
        let cT0D0 = self.cT0D0(eta);
        2.0 * self.cT1D0(eta) * self.cT0D1(eta) / cT0D0.powi(3) - self.cT1D1(eta) / cT0D0.powi(2)
    }
    fn c1T1D2(&self, eta: f64) -> f64 {
        let cT0D0 = self.cT0D0(eta);
        let cT0D1 = self.cT0D1(eta);
        let cT1D0 = self.cT1D0(eta);
        -6.0 * cT1D0 * cT0D1.powi(2) / cT0D0.powi(4)
            + 4.0 * self.cT1D1(eta) * cT0D1 / cT0D0.powi(3)
            + 2.0 * cT1D0 * self.cT0D2(eta) / cT0D0.powi(3)
            - self.cT1D2(eta) / cT0D0.powi(2)
    }
    fn c1T1D3(&self, eta: f64) -> f64 {
        let cT0D0 = self.cT0D0(eta);
        let cT0D1 = self.cT0D1(eta);
        let cT0D2 = self.cT0D2(eta);
        let cT1D0 = self.cT1D0(eta);
        let cT1D1 = self.cT1D1(eta);
        24.0 * cT1D0 * cT0D1.powi(3) / cT0D0.powi(5)
            - 18.0 * cT1D1 * cT0D1.powi(2) / cT0D0.powi(4)
            - 18.0 * cT1D0 * cT0D1 * cT0D2 / cT0D0.powi(4)
            + 6.0 * cT1D1 * cT0D2 / cT0D0.powi(3)
            + 6.0 * self.cT1D2(eta) * cT0D1 / cT0D0.powi(3)
            + 2.0 * cT1D0 * self.cT0D3(eta) / cT0D0.powi(3)
            - self.cT1D3(eta) / cT0D0.powi(2)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn cT0D0(&self, eta: f64) -> f64 {
        1.0 + 2.0 * self.m * (4.0 * eta - eta.powi(2)) / (1.0 - eta).powi(4)
            + (1.0 - self.m)
                * (20.0 * eta - 27.0 * eta.powi(2) + 12.0 * eta.powi(3) - 2.0 * eta.powi(4))
                / ((1.0 - eta) * (2.0 - eta)).powi(2)
    }
    fn cT0D1(&self, eta: f64) -> f64 {
        eta * 2.0
            * (2.0 * self.m * (-eta.powi(2) + 5.0 * eta + 2.0) / (1.0 - eta).powi(5)
                + (1.0 - self.m) * (eta.powi(3) + 6.0 * eta.powi(2) - 24.0 * eta + 20.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(3))
    }
    fn cT0D2(&self, eta: f64) -> f64 {
        eta.powi(2)
            * 6.0
            * (2.0 * self.m * (-eta.powi(2) + 6.0 * eta + 5.0) / (1.0 - eta).powi(6)
                + (1.0 - self.m)
                    * (-eta.powi(4) - 8.0 * eta.powi(3) + 48.0 * eta.powi(2) - 80.0 * eta + 44.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(4))
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
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn i1T0D0(&self, eta: f64) -> f64 {
        (A0[0] + self.m1 * A1[0] + self.m12 * A2[0])
            + (A0[1] + self.m1 * A1[1] + self.m12 * A2[1]) * eta
            + (A0[2] + self.m1 * A1[2] + self.m12 * A2[2]) * eta.powi(2)
            + (A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * eta.powi(3)
            + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * eta.powi(4)
            + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * eta.powi(5)
            + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * eta.powi(6)
    }
    fn i1T0D1(&self, eta: f64) -> f64 {
        (A0[1] + self.m1 * A1[1] + self.m12 * A2[1]) * eta
            + (A0[2] + self.m1 * A1[2] + self.m12 * A2[2]) * 2.0 * eta.powi(2)
            + (A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * 3.0 * eta.powi(3)
            + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 4.0 * eta.powi(4)
            + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 5.0 * eta.powi(5)
            + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 6.0 * eta.powi(6)
    }
    fn i1T0D2(&self, eta: f64) -> f64 {
        (A0[2] + self.m1 * A1[2] + self.m12 * A2[2]) * 2.0 * eta.powi(2)
            + (A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * 6.0 * eta.powi(3)
            + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 12.0 * eta.powi(4)
            + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 20.0 * eta.powi(5)
            + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 30.0 * eta.powi(6)
    }
    fn i1T0D3(&self, eta: f64) -> f64 {
        (A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * 6.0 * eta.powi(3)
            + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 24.0 * eta.powi(4)
            + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 60.0 * eta.powi(5)
            + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 120.0 * eta.powi(6)
    }
    fn i1T0D4(&self, eta: f64) -> f64 {
        (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 24.0 * eta.powi(4)
            + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 120.0 * eta.powi(5)
            + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 360.0 * eta.powi(6)
    }
    fn i1T1D0(&self, eta: f64) -> f64 {
        self.eta1
            * ((A0[1] + self.m1 * A1[1] + self.m12 * A2[1])
                + (A0[2] + self.m1 * A1[2] + self.m12 * A2[2]) * 2.0 * eta
                + (A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * 3.0 * eta.powi(2)
                + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 4.0 * eta.powi(3)
                + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 5.0 * eta.powi(4)
                + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 6.0 * eta.powi(5))
    }
    fn i1T1D1(&self, eta: f64) -> f64 {
        self.eta1
            * ((A0[1] + self.m1 * A1[1] + self.m12 * A2[1])
                + (A0[2] + self.m1 * A1[2] + self.m12 * A2[2]) * 4.0 * eta
                + (A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * 9.0 * eta.powi(2)
                + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 16.0 * eta.powi(3)
                + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 25.0 * eta.powi(4)
                + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 36.0 * eta.powi(5))
    }
    fn i1T1D2(&self, eta: f64) -> f64 {
        self.eta1
            * ((A0[2] + self.m1 * A1[2] + self.m12 * A2[2]) * 4.0 * eta
                + (A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * 18.0 * eta.powi(2)
                + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 48.0 * eta.powi(3)
                + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 100.0 * eta.powi(4)
                + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 180.0 * eta.powi(5))
    }
    fn i1T1D3(&self, eta: f64) -> f64 {
        self.eta1
            * ((A0[3] + self.m1 * A1[3] + self.m12 * A2[3]) * 18.0 * eta.powi(2)
                + (A0[4] + self.m1 * A1[4] + self.m12 * A2[4]) * 96.0 * eta.powi(3)
                + (A0[5] + self.m1 * A1[5] + self.m12 * A2[5]) * 300.0 * eta.powi(4)
                + (A0[6] + self.m1 * A1[6] + self.m12 * A2[6]) * 720.0 * eta.powi(5))
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn i2T0D0(&self, eta: f64) -> f64 {
        (B0[0] + self.m1 * B1[0] + self.m12 * B2[0])
            + (B0[1] + self.m1 * B1[1] + self.m12 * B2[1]) * eta
            + (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * eta.powi(2)
            + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * eta.powi(3)
            + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * eta.powi(4)
            + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * eta.powi(5)
            + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * eta.powi(6)
    }
    fn i2T0D1(&self, eta: f64) -> f64 {
        (B0[1] + self.m1 * B1[1] + self.m12 * B2[1]) * eta
            + (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * 2.0 * eta.powi(2)
            + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 3.0 * eta.powi(3)
            + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 4.0 * eta.powi(4)
            + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 5.0 * eta.powi(5)
            + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 6.0 * eta.powi(6)
    }
    fn i2T0D2(&self, eta: f64) -> f64 {
        (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * 2.0 * eta.powi(2)
            + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 6.0 * eta.powi(3)
            + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 12.0 * eta.powi(4)
            + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 20.0 * eta.powi(5)
            + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 30.0 * eta.powi(6)
    }
    fn i2T0D3(&self, eta: f64) -> f64 {
        (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 6.0 * eta.powi(3)
            + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 24.0 * eta.powi(4)
            + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 60.0 * eta.powi(5)
            + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 120.0 * eta.powi(6)
    }
    fn i2T0D4(&self, eta: f64) -> f64 {
        (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 24.0 * eta.powi(4)
            + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 120.0 * eta.powi(5)
            + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 360.0 * eta.powi(6)
    }
    fn i2T1D0(&self, eta: f64) -> f64 {
        self.eta1
            * ((B0[1] + self.m1 * B1[1] + self.m12 * B2[1])
                + (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * 2.0 * eta
                + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 3.0 * eta.powi(2)
                + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 4.0 * eta.powi(3)
                + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 5.0 * eta.powi(4)
                + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 6.0 * eta.powi(5))
    }
    fn i2T1D1(&self, eta: f64) -> f64 {
        self.eta1
            * ((B0[1] + self.m1 * B1[1] + self.m12 * B2[1])
                + (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * 4.0 * eta
                + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 9.0 * eta.powi(2)
                + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 16.0 * eta.powi(3)
                + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 25.0 * eta.powi(4)
                + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 36.0 * eta.powi(5))
    }
    fn i2T1D2(&self, eta: f64) -> f64 {
        self.eta1
            * ((B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * 4.0 * eta
                + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 18.0 * eta.powi(2)
                + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 48.0 * eta.powi(3)
                + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 100.0 * eta.powi(4)
                + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 180.0 * eta.powi(5))
    }
    fn i2T1D3(&self, eta: f64) -> f64 {
        self.eta1
            * ((B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 18.0 * eta.powi(2)
                + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 96.0 * eta.powi(3)
                + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 300.0 * eta.powi(4)
                + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 720.0 * eta.powi(5))
    }
}
const R: f64 = 8.314462618;
const NA: f64 = 6.02214076E23;
const PI: f64 = std::f64::consts::PI;
const A0: [f64; 7] = [
    0.9105631445,
    0.6361281449,
    2.6861347891,
    -26.547362491,
    97.759208784,
    -159.59154087,
    91.297774084,
];
const A1: [f64; 7] = [
    -0.3084016918,
    0.1860531159,
    -2.5030047259,
    21.419793629,
    -65.255885330,
    83.318680481,
    -33.746922930,
];
const A2: [f64; 7] = [
    -0.0906148351,
    0.4527842806,
    0.5962700728,
    -1.7241829131,
    -4.1302112531,
    13.776631870,
    -8.6728470368,
];
const B0: [f64; 7] = [
    0.7240946941,
    2.2382791861,
    -4.0025849485,
    -21.003576815,
    26.855641363,
    206.55133841,
    -355.60235612,
];
const B1: [f64; 7] = [
    -0.5755498075,
    0.6995095521,
    3.8925673390,
    -17.215471648,
    192.67226447,
    -161.82646165,
    -165.20769346,
];
const B2: [f64; 7] = [
    0.0976883116,
    -0.2557574982,
    -9.1558561530,
    20.642075974,
    -38.804430052,
    93.626774077,
    -29.666905585,
];
/// Romberg Numerical Differentiation
fn romberg_diff<F>(mut f: F, x: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    let mut fx: [[f64; 12]; 12] = [[0.0; 12]; 12];
    let mut h: f64 = 0.01 * x;
    fx[0][0] = (f(x + h) - f(x - h)) / 2.0 / h;
    for i in 1..12 {
        h /= 2.0;
        fx[0][i] = (f(x + h) - f(x - h)) / 2.0 / h;
        for j in 1..i + 1 {
            fx[j][i] = (fx[j - 1][i] * 4_f64.powi(j as i32) - fx[j - 1][i - 1])
                / (4_f64.powi(j as i32) - 1.0);
        }
        if (fx[i][i] / fx[i - 1][i - 1] - 1.0).abs() < 1E-10 {
            return fx[i][i];
        }
    }
    fx[11][11]
}
/// unit test
#[cfg(test)]
mod tests {
    #[test]
    fn test_pc_saft_pure() {}
}
