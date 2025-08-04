use super::{DispTerm, GiiPure, HsPure};
use super::{FRAC_NA_1E30, FRAC_RE30_NA, R};
use super::{PcSaftErr, PcSaftPure};
use crate::algorithms::{brent_zero, romberg_diff};
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_6};
#[cfg_attr(feature = "with_pyo3", pyclass)]
pub struct SPcSaftMix2 {
    x: [f64; 2],
    m: [f64; 2],
    m1: [f64; 2],
    sigma: [f64; 2],
    epsilon: [f64; 2],
    hs: HsPure,     // HsPure
    gii: GiiPure,   // GiiPure
    disp: DispTerm, // DispTerm
    m2e1s3_coef: [f64; 3],
    m2e2s3_coef: [f64; 3],
    // state
    temp: f64,
    rho_num: f64,
    eta0_coef: (f64, [f64; 2]),
    eta1_coef: (f64, [f64; 2]),
    eta2_coef: (f64, [f64; 2]),
    eta0_mu_k: (f64, [f64; 2]),
    is_single_phase: bool,
    // critical point
    omega1: Option<[f64; 2]>,
    temp_c: [f64; 2],
    pres_c: [f64; 2],
}
impl SPcSaftMix2 {
    pub fn new_fluid(
        x: [f64; 2],
        m: [f64; 2],
        sigma: [f64; 2],
        epsilon: [f64; 2],
        kij: f64,
    ) -> Self {
        let sigma3 = [sigma[0].powi(3), sigma[1].powi(3)];
        let m2e1s3_coef = [
            m[0].powi(2) * epsilon[0] * sigma3[0],
            m[1].powi(2) * epsilon[1] * sigma3[1],
            2.0 * (m[0] * m[1])
                * ((epsilon[0] * epsilon[1]).sqrt() * (1.0 - kij))
                * ((sigma[0] + sigma[1]).powi(3) / 8.0),
        ];
        let m2e2s3_coef = [
            (m[0] * epsilon[0]).powi(2) * sigma3[0],
            (m[1] * epsilon[1]).powi(2) * sigma3[1],
            2.0 * (m[0] * m[1])
                * (epsilon[0] * epsilon[1] * (1.0 - kij).powi(2))
                * ((sigma[0] + sigma[1]).powi(3) / 8.0),
        ];
        let m1 = [m[0] - 1.0, m[1] - 1.0];
        Self {
            x,
            m,
            m1,
            sigma,
            epsilon,
            hs: HsPure::new(x[0] * m[0] + x[1] * m[1]),
            gii: GiiPure::new(x[0] * m1[0] + x[1] * m1[1]),
            disp: DispTerm::new(
                x[0] * m[0] + x[1] * m[1],
                x[0].powi(2) * m2e1s3_coef[0]
                    + x[1].powi(2) * m2e1s3_coef[1]
                    + x[0] * x[1] * m2e1s3_coef[2],
                x[0].powi(2) * m2e2s3_coef[0]
                    + x[1].powi(2) * m2e2s3_coef[1]
                    + x[0] * x[1] * m2e2s3_coef[2],
            ),
            m2e1s3_coef,
            m2e2s3_coef,
            // state
            temp: 0.0,
            rho_num: 0.0,
            eta0_coef: (0.0, [0.0, 0.0]),
            eta1_coef: (0.0, [0.0, 0.0]),
            eta2_coef: (0.0, [0.0, 0.0]),
            eta0_mu_k: (0.0, [0.0, 0.0]),
            is_single_phase: true,
            // critical point
            omega1: None,
            temp_c: [0.0, 0.0],
            pres_c: [0.0, 0.0],
        }
    }
    fn new_fracs(&self, x: [f64; 2]) -> Self {
        Self {
            x,
            hs: HsPure::new(x[0] * self.m[0] + x[1] * self.m[1]),
            gii: GiiPure::new(x[0] * self.m1[0] + x[1] * self.m1[1]),
            disp: DispTerm::new(
                x[0] * self.m[0] + x[1] * self.m[1],
                x[0].powi(2) * self.m2e1s3_coef[0]
                    + x[1].powi(2) * self.m2e1s3_coef[1]
                    + x[0] * x[1] * self.m2e1s3_coef[2],
                x[0].powi(2) * self.m2e2s3_coef[0]
                    + x[1].powi(2) * self.m2e2s3_coef[1]
                    + x[0] * x[1] * self.m2e2s3_coef[2],
            ),
            ..*self
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
impl SPcSaftMix2 {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(x: [f64; 2], m: [f64; 2], sigma: [f64; 2], epsilon: [f64; 2], kij: f64) -> Self {
        Self::new_fluid(x, m, sigma, epsilon, kij)
    }
}
fn_tpz_flash_mix2!(SPcSaftMix2);
fn_tx_flash_mix2!(SPcSaftMix2);
fn_ty_flash_mix2!(SPcSaftMix2);
fn_tp_flash!(SPcSaftMix2);
fn_single_prop!(SPcSaftMix2);
impl SPcSaftMix2 {
    fn_calc_prop!();
    fn_check_derivatives!();
}
impl SPcSaftMix2 {
    fn eta0_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta0_coef.0 {
            self.eta0_coef = (
                temp,
                [
                    FRAC_PI_6
                        * self.m[0]
                        * (self.sigma[0] * (1.0 - 0.12 * (-3.0 * self.epsilon[0] / temp).exp()))
                            .powi(3),
                    FRAC_PI_6
                        * self.m[1]
                        * (self.sigma[1] * (1.0 - 0.12 * (-3.0 * self.epsilon[1] / temp).exp()))
                            .powi(3),
                ],
            );
        }
        self.x[0] * self.eta0_coef.1[0] + self.x[1] * self.eta0_coef.1[1]
    }
    fn eta1_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta1_coef.0 {
            let epsilon_temp_plus = [3.0 * self.epsilon[0] / temp, 3.0 * self.epsilon[1] / temp];
            self.eta1_coef = (
                temp,
                [
                    FRAC_PI_2
                        * self.m[0]
                        * (self.sigma[0] * (1.0 - 0.12 * (-epsilon_temp_plus[0]).exp())).powi(2)
                        * (self.sigma[0]
                            * (-0.12 * (-epsilon_temp_plus[0]).exp() * epsilon_temp_plus[0])),
                    FRAC_PI_2
                        * self.m[1]
                        * (self.sigma[1] * (1.0 - 0.12 * (-epsilon_temp_plus[1]).exp())).powi(2)
                        * (self.sigma[1]
                            * (-0.12 * (-epsilon_temp_plus[1]).exp() * epsilon_temp_plus[1])),
                ],
            );
        }
        self.x[0] * self.eta1_coef.1[0] + self.x[1] * self.eta1_coef.1[1]
    }
    fn eta2_coef(&mut self, temp: f64) -> f64 {
        if temp != self.eta2_coef.0 {
            let epsilon_temp_plus = [3.0 * self.epsilon[0] / temp, 3.0 * self.epsilon[1] / temp];
            self.eta2_coef = (
                temp,
                [
                    FRAC_PI_2
                        * self.m[0]
                        * self.sigma[0].powi(3)
                        * (2.0
                            * (1.0 - 0.12 * (-epsilon_temp_plus[0]).exp())
                            * (-0.12 * (-epsilon_temp_plus[0]).exp() * epsilon_temp_plus[0])
                                .powi(2)
                            + (1.0 - 0.12 * (-epsilon_temp_plus[0]).exp()).powi(2)
                                * (-0.12
                                    * (-epsilon_temp_plus[0]).exp()
                                    * epsilon_temp_plus[0]
                                    * (epsilon_temp_plus[0] - 2.0))),
                    FRAC_PI_2
                        * self.m[1]
                        * self.sigma[1].powi(3)
                        * (2.0
                            * (1.0 - 0.12 * (-epsilon_temp_plus[1]).exp())
                            * (-0.12 * (-epsilon_temp_plus[1]).exp() * epsilon_temp_plus[1])
                                .powi(2)
                            + (1.0 - 0.12 * (-epsilon_temp_plus[1]).exp()).powi(2)
                                * (-0.12
                                    * (-epsilon_temp_plus[1]).exp()
                                    * epsilon_temp_plus[1]
                                    * (epsilon_temp_plus[1] - 2.0))),
                ],
            );
        }
        self.x[0] * self.eta2_coef.1[0] + self.x[1] * self.eta2_coef.1[1]
    }
    fn eta0_mu_k(&mut self, temp: f64) -> [f64; 2] {
        if temp != self.eta0_mu_k.0 {
            self.eta0_mu_k = (
                temp,
                [
                    FRAC_PI_6
                        * self.m[0]
                        * (self.sigma[0] * (1.0 - 0.12 * (-3.0 * self.epsilon[0] / temp).exp()))
                            .powi(3),
                    FRAC_PI_6
                        * self.m[1]
                        * (self.sigma[1] * (1.0 - 0.12 * (-3.0 * self.epsilon[1] / temp).exp()))
                            .powi(3),
                ],
            );
        }
        self.eta0_mu_k.1
    }
    fn r_t0d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d0(eta) + self.gii.lngii_t0d0(eta) + self.disp.t0d0(temp, rho_num, eta)
    }
    fn r_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d1(eta) + self.gii.lngii_t0d1(eta) + self.disp.t0d1(temp, rho_num, eta)
    }
    fn r_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d2(eta) + self.gii.lngii_t0d2(eta) + self.disp.t0d2(temp, rho_num, eta)
    }
    fn r_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d3(eta) + self.gii.lngii_t0d3(eta) + self.disp.t0d3(temp, rho_num, eta)
    }
    fn r_t0d4(&mut self, temp: f64, rho_num: f64) -> f64 {
        let eta = self.eta0_coef(temp) * rho_num;
        self.hs.t0d4(eta) + self.gii.lngii_t0d4(eta) + self.disp.t0d4(temp, rho_num, eta)
    }
    fn r_t1d0(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d0(eta, eta1)
            + self.gii.lngii_t1d0(eta, eta1)
            + self.disp.t1d0(temp, rho_num, eta, eta1)
    }
    fn r_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d1(eta, eta1)
            + self.gii.lngii_t1d1(eta, eta1)
            + self.disp.t1d1(temp, rho_num, eta, eta1)
    }
    fn r_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d2(eta, eta1)
            + self.gii.lngii_t1d2(eta, eta1)
            + self.disp.t1d2(temp, rho_num, eta, eta1)
    }
    fn r_t1d3(&mut self, temp: f64, rho_num: f64) -> f64 {
        let (eta, eta1) = (
            self.eta0_coef(temp) * rho_num,
            self.eta1_coef(temp) * rho_num,
        );
        self.hs.t1d3(eta, eta1)
            + self.gii.lngii_t1d3(eta, eta1)
            + self.disp.t1d3(temp, rho_num, eta, eta1)
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
    }
    fn ln_phi(&mut self, temp: f64, pres: f64, eta0_guess: f64) -> [f64; 2] {
        let eta0_coef = self.eta0_coef(temp);
        let rho_num = self.calc_density(temp, pres, eta0_guess / eta0_coef);
        let eta = eta0_coef * rho_num;
        let eta_k = self.eta0_mu_k(temp);
        let m2e1s3_k = [
            2.0 * self.x[0] * self.m2e1s3_coef[0] + self.x[1] * self.m2e1s3_coef[2],
            2.0 * self.x[1] * self.m2e1s3_coef[1] + self.x[0] * self.m2e1s3_coef[2],
        ];
        let m2e2s3_k = [
            2.0 * self.x[0] * self.m2e2s3_coef[0] + self.x[1] * self.m2e2s3_coef[2],
            2.0 * self.x[1] * self.m2e2s3_coef[1] + self.x[0] * self.m2e2s3_coef[2],
        ];
        let t0d1 = self.r_t0d1(temp, rho_num);
        let ln_phi: Vec<f64> = self
            .hs
            .mu_k(eta, rho_num, &self.m, &eta_k)
            .zip(self.gii.lngii_mu_k(eta, rho_num, &self.m1, &eta_k))
            .zip(
                self.disp
                    .mu_k(temp, rho_num, eta, &self.m, &eta_k, &m2e1s3_k, &m2e2s3_k),
            )
            .map(|((hs, gii), disp)| hs + gii + disp - (1.0 + t0d1).ln())
            .collect();
        [ln_phi[0], ln_phi[1]]
    }
    fn calc_ln_k(&mut self, temp: f64, pres: f64, x: [f64; 2], y: [f64; 2]) -> [f64; 2] {
        let ln_phi_l = self.new_fracs(x).ln_phi(temp, pres, 0.5);
        let ln_phi_v = self.new_fracs(y).ln_phi(temp, pres, 1e-10);
        [ln_phi_l[0] - ln_phi_v[0], ln_phi_l[1] - ln_phi_v[1]]
    }
    fn guess_ps(&mut self, temp: f64) -> [f64; 2] {
        if self.omega1.is_none() {
            let mut fluid = [
                PcSaftPure::new_fluid(self.m[0], self.sigma[0], self.epsilon[0]),
                PcSaftPure::new_fluid(self.m[1], self.sigma[1], self.epsilon[1]),
            ];
            fluid[0].c_flash().unwrap();
            fluid[1].c_flash().unwrap();
            self.temp_c = [fluid[0].T().unwrap(), fluid[1].T().unwrap()];
            self.pres_c = [fluid[0].p().unwrap(), fluid[1].p().unwrap()];
            fluid[0].t_flash(0.7 * self.temp_c[0]).unwrap();
            fluid[1].t_flash(0.7 * self.temp_c[1]).unwrap();
            self.omega1 = Some([
                -(fluid[0].p_s().unwrap() / self.pres_c[0]).log10(),
                -(fluid[1].p_s().unwrap() / self.pres_c[1]).log10(),
            ]);
        }
        [
            self.pres_c[0]
                * 10_f64.powf(7.0 / 3.0 * self.omega1.unwrap()[0] * (1.0 - self.temp_c[0] / temp)),
            self.pres_c[1]
                * 10_f64.powf(7.0 / 3.0 * self.omega1.unwrap()[1] * (1.0 - self.temp_c[1] / temp)),
        ]
    }
}
#[cfg(test)]
mod tests {
    #[test]
    fn test_s_pc_saft_mix2() {
        let mut fluids = super::SPcSaftMix2::new_fluid(
            [0.5, 0.5],
            [1.0000, 1.6069],
            [3.7039, 3.5206],
            [150.03, 191.42],
            0.0,
        );
        let temp = 199.92;
        // test tpz_flash
        let p_data = [
            3.62e5, 4.42e5, 6.8e5, 10.9e5, 17e5, 23.8e5, 34e5, 40.8e5, 47.65e5, 48.9e5, 49.4e5,
            49.8e5, 50.35e5,
        ];
        [
            (0.0291, 0.3962),
            (0.0452, 0.5041),
            (0.0933, 0.6756),
            (0.1767, 0.7959),
            (0.3018, 0.8682),
            (0.4419, 0.9058),
            (0.6475, 0.9357),
            (0.7755, 0.9488),
            (0.8909, 0.9597),
            (0.9103, 0.9614),
            (0.9179, 0.9620),
            (0.9240, 0.9624),
            (0.9324, 0.9628),
        ]
        .iter()
        .zip(p_data.iter())
        .map(|(&(x, y), &p)| {
            let xy = fluids.tpz_flash(temp, p).unwrap();
            assert_eq!((xy[0] * 1e4).round() / 1e4, x);
            assert_eq!((xy[1] * 1e4).round() / 1e4, y);
        })
        .count();
        // test tx_flash
        let x_data = [
            0.0214, 0.0512, 0.1039, 0.1875, 0.31, 0.4526, 0.6601, 0.7852, 0.8942, 0.9126, 0.9175,
            0.9222, 0.9319,
        ];
        [
            (3.24e5, 0.3257),
            (4.72e5, 0.5350),
            (7.33e5, 0.6985),
            (11.43e5, 0.8052),
            (17.40e5, 0.8711),
            (24.32e5, 0.9079),
            (34.64e5, 0.9371),
            (41.34e5, 0.9498),
            (47.86e5, 0.9600),
            (49.05e5, 0.9616),
            (49.37e5, 0.9619),
            (49.68e5, 0.9622),
            (50.32e5, 0.9627),
        ]
        .iter()
        .zip(x_data.iter())
        .map(|(&(p, y), &x)| {
            let py = fluids.tx_flash(temp, x).unwrap();
            assert_eq!((py[0] / 1e3).round() * 1e3, p);
            assert_eq!((py[1] * 1e4).round() / 1e4, y);
        })
        .count();
        // test ty_flash
        let y_data = [
            0.3005, 0.5098, 0.68, 0.7957, 0.8679, 0.9052, 0.9337, 0.9461, 0.9562, 0.9584, 0.9578,
            0.9575, 0.9577,
        ];
        [
            (3.12e5, 0.0190),
            (4.47e5, 0.0462),
            (6.90e5, 0.0952),
            (10.89e5, 0.1764),
            (16.97e5, 0.3011),
            (23.65e5, 0.4388),
            (33.10e5, 0.6299),
            (39.25e5, 0.7474),
        ]
        .iter()
        .zip(y_data.iter())
        .map(|(&(p, x), &y)| {
            let px = fluids.ty_flash(temp, y).unwrap();
            assert_eq!((px[0] / 1e3).round() * 1e3, p);
            assert_eq!((px[1] * 1e4).round() / 1e4, x);
        })
        .count();
    }
}
