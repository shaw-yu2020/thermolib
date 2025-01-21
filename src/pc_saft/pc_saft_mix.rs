use super::{CTerm, I1Term, I2Term};
use crate::algorithms::romberg_diff;
use crate::f64consts::{FRAC_NA_1E30, FRAC_PI_6, PI, R};
#[allow(non_snake_case)]
pub struct PcSaftMix {
    params: Vec<PcSaftParams>,
    m: f64,
    m2e1s3: f64,
    m2e2s3: f64,
    sum_xmd1: (f64, f64),
    sum_xmd2: (f64, f64),
    sum_xmd3: (f64, f64),
    i1: I1Term, // I1Term
    i2: I2Term, // I2Term
    c: CTerm,   // Cterm
}
#[allow(non_snake_case)]
impl PcSaftMix {
    pub fn new_fluid(x: Vec<f64>, m: Vec<f64>, sigma: Vec<f64>, epsilon: Vec<f64>) -> Self {
        let mut params = Vec::with_capacity(x.len());
        if x.len() == m.len() && x.len() == sigma.len() && x.len() == epsilon.len() {
            for i in 0..x.len() {
                params.push(PcSaftParams {
                    x: x[i],
                    m: m[i],
                    sigma: sigma[i],
                    epsillon: epsilon[i],
                    d: (0.0, 0.0),
                });
            }
        }
        let m = params.iter().map(|param| param.xm()).sum();
        Self {
            m2e1s3: params
                .iter()
                .map(|i| {
                    params
                        .iter()
                        .map(|j| {
                            (i.xm() * j.xm())
                                * (i.epsillon * j.epsillon).sqrt()
                                * (i.sigma + j.sigma).powi(3)
                        })
                        .sum::<f64>()
                })
                .sum::<f64>()
                / 8.0,
            m2e2s3: params
                .iter()
                .map(|i| {
                    params
                        .iter()
                        .map(|j| {
                            i.xm() * j.xm() * i.epsillon * j.epsillon * (i.sigma + j.sigma).powi(3)
                        })
                        .sum::<f64>()
                })
                .sum::<f64>()
                / 8.0,
            params,
            m,
            sum_xmd1: (0.0, 0.0),
            sum_xmd2: (0.0, 0.0),
            sum_xmd3: (0.0, 0.0),
            i1: I1Term::new(m), // I1Term
            i2: I2Term::new(m), // I1Term
            c: CTerm::new(m),   // CTerm
        }
    }
    pub fn calc_density(&mut self, temp: f64, p: f64, eta_guess: f64) -> f64 {
        let mut rho_num = eta_guess / (FRAC_PI_6 * self.sum_xmd3(temp)); // rho_num_guess
        let (mut p_diff, mut rho_num_diff, mut val_Dp_Drho_T);
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
impl PcSaftMix {
    fn calc_p(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp) * rho_num * (1.0 + self.calc_rT0D1(temp, rho_num))
    }
    fn calc_Dp_Drho_T(&mut self, temp: f64, rho_num: f64) -> f64 {
        (R / FRAC_NA_1E30 * temp)
            * (1.0 + 2.0 * self.calc_rT0D1(temp, rho_num) + self.calc_rT0D2(temp, rho_num))
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn calc_rT0D0(&mut self, temp: f64, dens: f64) -> f64 {
        self.hsT0D0(temp, dens) - self.lngiiT0D0(temp, dens) + self.dispT0D0(temp, dens)
    }
    fn calc_rT0D1(&mut self, temp: f64, dens: f64) -> f64 {
        self.hsT0D1(temp, dens) - self.lngiiT0D1(temp, dens) + self.dispT0D1(temp, dens)
    }
    fn calc_rT0D2(&mut self, temp: f64, dens: f64) -> f64 {
        self.hsT0D2(temp, dens) - self.lngiiT0D2(temp, dens) + self.dispT0D2(temp, dens)
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn hsT0D0(&mut self, temp: f64, dens: f64) -> f64 {
        let dens_plus = FRAC_PI_6 * dens;
        let zeta3_plus = 1.0 - dens_plus * self.sum_xmd3(temp);
        3.0 * dens_plus * self.sum_xmd1(temp) * self.sum_xmd2(temp) / zeta3_plus
            + self.sum_xmd2(temp).powi(3) / self.sum_xmd3(temp)
                * (dens_plus / zeta3_plus.powi(2) + zeta3_plus.ln() / self.sum_xmd3(temp))
            - self.m * zeta3_plus.ln()
    }
    fn hsT0D1(&mut self, temp: f64, dens: f64) -> f64 {
        let dens_plus = FRAC_PI_6 * dens;
        let zeta3 = dens_plus * self.sum_xmd3(temp);
        3.0 * dens_plus * self.sum_xmd1(temp) * self.sum_xmd2(temp) / (1.0 - zeta3).powi(2)
            + dens_plus.powi(2) * (self.sum_xmd2(temp) / (1.0 - zeta3)).powi(3) * (3.0 - zeta3)
            + self.m * zeta3 / (1.0 - zeta3)
    }
    fn hsT0D2(&mut self, temp: f64, dens: f64) -> f64 {
        let dens_plus = FRAC_PI_6 * dens;
        let zeta3 = dens_plus * self.sum_xmd3(temp);
        dens_plus.powi(2) / (1.0 - zeta3).powi(3)
            * (6.0 * self.sum_xmd1(temp) * self.sum_xmd2(temp) * self.sum_xmd3(temp)
                + self.sum_xmd2(temp).powi(3) / (1.0 - zeta3) * (3.0 + 4.0 * zeta3 - zeta3.powi(2)))
            + self.m * (zeta3 / (1.0 - zeta3)).powi(2)
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn lngiiT0D0(&mut self, temp: f64, dens: f64) -> f64 {
        let dens_plus = FRAC_PI_6 * dens;
        let giiT0D0d0 = (1.0 - dens_plus * self.sum_xmd3(temp)).recip();
        let zeta2plus = dens_plus * self.sum_xmd2(temp) * giiT0D0d0;
        let giiT0D0d1 = giiT0D0d0 * 1.5 * zeta2plus;
        let giiT0D0d2 = giiT0D0d0 * 0.5 * zeta2plus.powi(2);
        self.params
            .iter_mut()
            .map(|params| {
                params.xm1()
                    * (giiT0D0d0 + giiT0D0d1 * params.d(temp) + giiT0D0d2 * params.d(temp).powi(2))
                        .ln()
            })
            .sum()
    }
    fn lngiiT0D1(&mut self, temp: f64, dens: f64) -> f64 {
        let dens_plus = FRAC_PI_6 * dens;
        let giiT0D0d0 = (1.0 - dens_plus * self.sum_xmd3(temp)).recip();
        let zeta2plus = dens_plus * self.sum_xmd2(temp) * giiT0D0d0;
        let zeta3plus = dens_plus * self.sum_xmd3(temp) * giiT0D0d0;
        let giiT0D0d1 = giiT0D0d0 * 1.5 * zeta2plus;
        let giiT0D0d2 = giiT0D0d0 * 0.5 * zeta2plus.powi(2);
        let giiT0D1d0 = giiT0D0d0 * zeta3plus;
        let giiT0D1d1 = giiT0D0d1 * (1.0 + 2.0 * zeta3plus);
        let giiT0D1d2 = giiT0D0d2 * (2.0 + 3.0 * zeta3plus);
        self.params
            .iter_mut()
            .map(|params| {
                params.xm1()
                    * (giiT0D1d0 + giiT0D1d1 * params.d(temp) + giiT0D1d2 * params.d(temp).powi(2))
                    / (giiT0D0d0 + giiT0D0d1 * params.d(temp) + giiT0D0d2 * params.d(temp).powi(2))
            })
            .sum()
    }
    fn lngiiT0D2(&mut self, temp: f64, dens: f64) -> f64 {
        let dens_plus = FRAC_PI_6 * dens;
        let giiT0D0d0 = (1.0 - dens_plus * self.sum_xmd3(temp)).recip();
        let zeta2plus = dens_plus * self.sum_xmd2(temp) * giiT0D0d0;
        let zeta3plus = dens_plus * self.sum_xmd3(temp) * giiT0D0d0;
        let giiT0D0d1 = giiT0D0d0 * 1.5 * zeta2plus;
        let giiT0D0d2 = giiT0D0d0 * 0.5 * zeta2plus.powi(2);
        let giiT0D1d0 = giiT0D0d0 * zeta3plus;
        let giiT0D1d1 = giiT0D0d1 * (1.0 + 2.0 * zeta3plus);
        let giiT0D1d2 = giiT0D0d2 * (2.0 + 3.0 * zeta3plus);
        let giiT0D2d0 = giiT0D0d0 * 2.0 * zeta3plus.powi(2);
        let giiT0D2d1 = giiT0D0d1 * zeta3plus * (4.0 + 6.0 * zeta3plus);
        let giiT0D2d2 = giiT0D0d2 * (2.0 + 12.0 * zeta3plus * (1.0 + zeta3plus));
        self.params
            .iter_mut()
            .map(|params| {
                let d2 = params.d(temp).powi(2);
                params.xm1() * (giiT0D2d0 + giiT0D2d1 * params.d(temp) + giiT0D2d2 * d2)
                    / (giiT0D0d0 + giiT0D0d1 * params.d(temp) + giiT0D0d2 * d2)
                    - params.xm1()
                        * ((giiT0D1d0 + giiT0D1d1 * params.d(temp) + giiT0D1d2 * d2)
                            / (giiT0D0d0 + giiT0D0d1 * params.d(temp) + giiT0D0d2 * d2))
                            .powi(2)
            })
            .sum()
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn dispT0D0(&mut self, temp: f64, dens: f64) -> f64 {
        let eta = FRAC_PI_6 * dens * self.sum_xmd3(temp);
        -PI * dens
            * (2.0 * self.i1T0D0(eta) * self.m2e1s3 / temp
                + self.m * self.c1T0D0(eta) * self.i2T0D0(eta) * self.m2e2s3 / temp.powi(2))
    }
    fn dispT0D1(&mut self, temp: f64, dens: f64) -> f64 {
        let eta = FRAC_PI_6 * dens * self.sum_xmd3(temp);
        -PI * dens
            * (2.0 * self.m2e1s3 / temp * (self.i1T0D0(eta) + self.i1T0D1(eta))
                + self.m * self.m2e2s3 / temp.powi(2)
                    * (self.c1T0D0(eta) * self.i2T0D0(eta)
                        + self.c1T0D1(eta) * self.i2T0D0(eta)
                        + self.c1T0D0(eta) * self.i2T0D1(eta)))
    }
    fn dispT0D2(&mut self, temp: f64, dens: f64) -> f64 {
        let eta = FRAC_PI_6 * dens * self.sum_xmd3(temp);
        -PI * dens
            * (2.0 * self.m2e1s3 / temp * (2.0 * self.i1T0D1(eta) + self.i1T0D2(eta))
                + self.m * self.m2e2s3 / temp.powi(2)
                    * (2.0 * self.c1T0D1(eta) * self.i2T0D0(eta)
                        + 2.0 * self.c1T0D0(eta) * self.i2T0D1(eta)
                        + self.c1T0D2(eta) * self.i2T0D0(eta)
                        + 2.0 * self.c1T0D1(eta) * self.i2T0D1(eta)
                        + self.c1T0D0(eta) * self.i2T0D2(eta)))
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn c1T0D0(&mut self, eta: f64) -> f64 {
        self.cT0D0(eta).recip()
    }
    fn c1T0D1(&mut self, eta: f64) -> f64 {
        -self.cT0D1(eta) / self.cT0D0(eta).powi(2)
    }
    fn c1T0D2(&mut self, eta: f64) -> f64 {
        2.0 * self.cT0D1(eta).powi(2) / self.cT0D0(eta).powi(3)
            - self.cT0D2(eta) / self.cT0D0(eta).powi(2)
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn cT0D0(&mut self, eta: f64) -> f64 {
        self.c.eta0(eta)
    }
    fn cT0D1(&mut self, eta: f64) -> f64 {
        eta * self.c.eta1(eta)
    }
    fn cT0D2(&mut self, eta: f64) -> f64 {
        eta.powi(2) * self.c.eta2(eta)
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn i1T0D0(&mut self, eta: f64) -> f64 {
        self.i1.eta0(eta)
    }
    fn i1T0D1(&mut self, eta: f64) -> f64 {
        eta * self.i1.eta1(eta)
    }
    fn i1T0D2(&mut self, eta: f64) -> f64 {
        eta.powi(2) * self.i1.eta2(eta)
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
    fn i2T0D0(&mut self, eta: f64) -> f64 {
        self.i2.eta0(eta)
    }
    fn i2T0D1(&mut self, eta: f64) -> f64 {
        eta * self.i2.eta1(eta)
    }
    fn i2T0D2(&mut self, eta: f64) -> f64 {
        eta.powi(2) * self.i2.eta2(eta)
    }
}
impl PcSaftMix {
    fn sum_xmd1(&mut self, temp: f64) -> f64 {
        if temp != self.sum_xmd1.0 {
            self.sum_xmd1 = (
                temp,
                self.params
                    .iter_mut()
                    .map(|params| params.xm() * params.d(temp))
                    .sum::<f64>(),
            )
        }
        self.sum_xmd1.1
    }
    fn sum_xmd2(&mut self, temp: f64) -> f64 {
        if temp != self.sum_xmd2.0 {
            self.sum_xmd2 = (
                temp,
                self.params
                    .iter_mut()
                    .map(|params| params.xm() * params.d(temp).powi(2))
                    .sum::<f64>(),
            )
        }
        self.sum_xmd2.1
    }
    fn sum_xmd3(&mut self, temp: f64) -> f64 {
        if temp != self.sum_xmd3.0 {
            self.sum_xmd3 = (
                temp,
                self.params
                    .iter_mut()
                    .map(|params| params.xm() * params.d(temp).powi(3))
                    .sum::<f64>(),
            )
        }
        self.sum_xmd3.1
    }
}
struct PcSaftParams {
    pub x: f64,
    pub m: f64,
    pub sigma: f64,
    pub epsillon: f64,
    /// Cached Variables
    d: (f64, f64), // (temp,value)
}
impl PcSaftParams {
    fn xm(&self) -> f64 {
        self.x * self.m
    }
    fn xm1(&self) -> f64 {
        self.x * (self.m - 1.0)
    }
    fn d(&mut self, temp: f64) -> f64 {
        if temp != self.d.0 {
            self.d = (
                temp,
                self.sigma * (1.0 - 0.12 * (-3.0 * self.epsillon / temp).exp()),
            );
        }
        self.d.1
    }
}
impl PcSaftMix {
    pub fn check_derivatives(&mut self, print_val: bool, temp: f64, dens: f64) {
        if print_val {
            println!(
                "[rT0D0 == rT0D0] calc_rT0D0() ={}",
                self.calc_rT0D0(temp, dens)
            );
        }
        let compare_val = |val_calc: f64, val_diff: f64| {
            assert_eq!(
                &val_calc.abs().to_string()[0..12],
                &val_diff.abs().to_string()[0..12]
            )
        };
        // derivative for density
        let val_calc = self.calc_rT0D1(temp, dens) / dens;
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT0D0(temp, rhox), dens);
        if print_val {
            println!("[rT0D1 == rT0D1] calc_rT0D1() ={}", val_calc);
            println!("[rT0D0 -> rT0D1] romberg_diff ={}", val_diff);
        } else {
            compare_val(val_calc, val_diff);
        }
        // derivative for density+density
        let val_calc = self.calc_rT0D2(temp, dens) / dens.powi(2);
        let val_diff = romberg_diff(|rhox: f64| self.calc_rT0D1(temp, rhox) / rhox, dens);
        if print_val {
            println!("[rT0D2 == rT0D2] calc_rT0D2() ={}", val_calc);
            println!("[rT0D1 -> rT0D2] romberg_diff ={}", val_diff);
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
    fn test_pc_saft_mix() {
        let (m_c11, sigma_c11, epsilon_c11) = (4.9082, 3.8893, 248.82);
        let (m_c12, sigma_c12, epsilon_c12) = (5.3060, 3.8959, 249.21);
        let mut c11_c12 = PcSaftMix::new_fluid(
            vec![0.5, 0.5],
            vec![m_c11, m_c12],
            vec![sigma_c11, sigma_c12],
            vec![epsilon_c11, epsilon_c12],
        );
        c11_c12.check_derivatives(false, 300.0, 0.003);
    }
}
