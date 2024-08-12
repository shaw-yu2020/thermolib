use thiserror::Error;
#[derive(Debug, Error)]
enum PcSaftPureErr {
    #[error("tp_flash diverge")]
    NotConvForTP,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
    #[error("property not in single phase")]
    NotInSinglePhase,
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
    eta: f64,
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
            T: 0.0,
            rho_num: 0.0,
            rhov_num: 0.0,
            rhol_num: 0.0,
            is_single_phase: true,
            m1: (m - 1.0) / m,
            m12: (m - 1.0) / m * (m - 2.0) / m,
            eta: 0.0,
            m2e1s3: 0.0,
            m2e2s3: 0.0,
        }
    }
    pub fn td_unchecked(&mut self, T: f64, rho_mol: f64) {
        self.set_temperature_and_number_density(T, rho_mol * NA / 1E30);
        self.is_single_phase = true;
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
        self.set_temperature_and_number_density(T, rho_num);
        (self.rho_num * 1E30 / NA) * R * T * (1.0 + self.rT0D1())
    }
    fn calc_Dp_Drho_T(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        1E30 / NA * R * T * (1.0 + 2.0 * self.rT0D1() + self.rT0D2())
    }
    fn calc_lnphi(&mut self, T: f64, rho_num: f64) -> f64 {
        self.set_temperature_and_number_density(T, rho_num);
        self.rT0D0() + self.rT0D1() - (1.0 + self.rT0D1()).ln()
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
            self.eta = PI / 6.0 * rho_num * self.m * d.powi(3);
        } else if rho_num != self.rho_num {
            self.rho_num = rho_num;
            let d = self.sigma * (1.0 - 0.12 * (-3.0 * self.epsilon / T).exp());
            self.eta = PI / 6.0 * rho_num * self.m * d.powi(3);
        }
    }
    fn rT0D0(&self) -> f64 {
        self.hcT0D0(self.eta) + self.rho_num * self.dispT0D0(self.eta)
    }
    fn rT0D1(&self) -> f64 {
        self.hcT0D1(self.eta) + self.rho_num * self.dispT0D1(self.eta)
    }
    fn rT0D2(&self) -> f64 {
        self.hcT0D2(self.eta) + self.rho_num * self.dispT0D2(self.eta)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn hcT0D0(&self, eta: f64) -> f64 {
        self.m * self.hsT0D0(eta) + (1.0 - self.m) * self.giiT0D0(eta).ln()
    }
    fn hcT0D1(&self, eta: f64) -> f64 {
        self.m * self.hsT0D1(eta) + (1.0 - self.m) / self.giiT0D0(eta) * self.giiT0D1(eta)
    }
    fn hcT0D2(&self, eta: f64) -> f64 {
        let giiT0D0: f64 = self.giiT0D0(eta);
        self.m * self.hsT0D2(eta)
            + (1.0 - self.m)
                * (-1.0 / giiT0D0.powi(2) * self.giiT0D1(eta).powi(2)
                    + 1.0 / giiT0D0 * self.giiT0D2(eta))
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn hsT0D0(&self, eta: f64) -> f64 {
        eta * (4.0 - 3.0 * eta) / (1.0 - eta).powi(2)
    }
    fn hsT0D1(&self, eta: f64) -> f64 {
        eta * (4.0 - 2.0 * eta) / (1.0 - eta).powi(3)
    }
    fn hsT0D2(&self, eta: f64) -> f64 {
        eta.powi(2) * (10.0 - 4.0 * eta) / (1.0 - eta).powi(4)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn giiT0D0(&self, eta: f64) -> f64 {
        (1.0 - 0.5 * eta) / (1.0 - eta).powi(3)
    }
    fn giiT0D1(&self, eta: f64) -> f64 {
        eta * (2.5 - eta) / (1.0 - eta).powi(4)
    }
    fn giiT0D2(&self, eta: f64) -> f64 {
        eta.powi(2) * (9.0 - 3.0 * eta) / (1.0 - eta).powi(5)
    }
}
#[allow(non_snake_case)]
impl PcSaftPure {
    fn dispT0D0(&self, eta: f64) -> f64 {
        -2.0 * PI * self.i1T0D0(eta) * self.m2e1s3
            - PI * self.m * self.c1T0D0(eta) * self.i2T0D0(eta) * self.m2e2s3
    }
    fn dispT0D1(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let i2T0D0 = self.i2T0D0(eta);
        -2.0 * PI * self.m2e1s3 * (self.i1T0D0(eta) + self.i1T0D1(eta))
            - PI * self.m
                * self.m2e2s3
                * (c1T0D0 * i2T0D0 + self.c1T0D1(eta) * i2T0D0 + c1T0D0 * self.i2T0D1(eta))
    }
    fn dispT0D2(&self, eta: f64) -> f64 {
        let c1T0D0 = self.c1T0D0(eta);
        let c1T0D1 = self.c1T0D1(eta);
        let i2T0D0 = self.i2T0D0(eta);
        let i2T0D1 = self.i2T0D1(eta);
        -2.0 * PI * self.m2e1s3 * (2.0 * self.i1T0D1(eta) + self.i1T0D2(eta))
            - PI * self.m
                * self.m2e2s3
                * (2.0 * c1T0D1 * i2T0D0
                    + 2.0 * c1T0D0 * i2T0D1
                    + self.c1T0D2(eta) * i2T0D0
                    + 2.0 * c1T0D1 * i2T0D1
                    + c1T0D0 * self.i2T0D2(eta))
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
        eta * (4.0 * self.m * (-eta.powi(2) + 5.0 * eta + 2.0) / (1.0 - eta).powi(5)
            + 2.0 * (1.0 - self.m) * (eta.powi(3) + 6.0 * eta.powi(2) - 24.0 * eta + 20.0)
                / ((1.0 - eta) * (2.0 - eta)).powi(3))
    }
    fn cT0D2(&self, eta: f64) -> f64 {
        eta.powi(2)
            * (12.0 * self.m * (-eta.powi(2) + 6.0 * eta + 5.0) / (1.0 - eta).powi(6)
                + (6.0 * (1.0 - self.m))
                    * (-eta.powi(4) - 8.0 * eta.powi(3) + 48.0 * eta.powi(2) - 80.0 * eta + 44.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(4))
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
impl PcSaftPure {}
