use crate::algorithms::romberg_diff;
use std::f64::consts::{FRAC_PI_6, PI};
const R: f64 = 8.314462618;
const FRAC_NA_1E30: f64 = 6.02214076E-7; // const NA: f64 = 6.02214076E23;
#[allow(non_snake_case)]
pub struct PcSaftMix {
    params: Vec<PcSaftParams>,
    m: f64,
    m1: f64,  // (m-1)/m
    m12: f64, // (m-1)/m * (m-2)/m
    m2e1s3: f64,
    m2e2s3: f64,
    sum_xmd1: (f64, f64),
    sum_xmd2: (f64, f64),
    sum_xmd3: (f64, f64),
    cT0D0: (f64, f64),
    cT0D1: (f64, f64),
    i2T0D0: (f64, f64),
    i2T0D1: (f64, f64),
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
            m1: (m - 1.0) / m,                      // (m-1)/m
            m12: (m - 1.0) * (m - 2.0) / m.powi(2), // (m-1)/m * (m-2)/m
            sum_xmd1: (0.0, 0.0),
            sum_xmd2: (0.0, 0.0),
            sum_xmd3: (0.0, 0.0),
            cT0D0: (0.0, 0.0),
            cT0D1: (0.0, 0.0),
            i2T0D0: (0.0, 0.0),
            i2T0D1: (0.0, 0.0),
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
        if eta != self.cT0D0.0 {
            self.cT0D0 = (
                eta,
                1.0 + 2.0 * self.m * (4.0 * eta - eta.powi(2)) / (1.0 - eta).powi(4)
                    + (1.0 - self.m)
                        * (20.0 * eta - 27.0 * eta.powi(2) + 12.0 * eta.powi(3)
                            - 2.0 * eta.powi(4))
                        / ((1.0 - eta) * (2.0 - eta)).powi(2),
            )
        }
        self.cT0D0.1
    }
    fn cT0D1(&mut self, eta: f64) -> f64 {
        if eta != self.cT0D1.0 {
            self.cT0D1 = (
                eta,
                (eta * 2.0)
                    * (2.0 * self.m * (-eta.powi(2) + 5.0 * eta + 2.0) / (1.0 - eta).powi(5)
                        + (1.0 - self.m) * (eta.powi(3) + 6.0 * eta.powi(2) - 24.0 * eta + 20.0)
                            / ((1.0 - eta) * (2.0 - eta)).powi(3)),
            )
        }
        self.cT0D1.1
    }
    fn cT0D2(&self, eta: f64) -> f64 {
        (eta.powi(2) * 6.0)
            * (2.0 * self.m * (-eta.powi(2) + 6.0 * eta + 5.0) / (1.0 - eta).powi(6)
                + (1.0 - self.m)
                    * (-eta.powi(4) - 8.0 * eta.powi(3) + 48.0 * eta.powi(2) - 80.0 * eta + 44.0)
                    / ((1.0 - eta) * (2.0 - eta)).powi(4))
    }
}
#[allow(non_snake_case)]
impl PcSaftMix {
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
#[allow(non_snake_case)]
impl PcSaftMix {
    fn i2T0D0(&mut self, eta: f64) -> f64 {
        if eta != self.i2T0D0.0 {
            self.i2T0D0 = (
                eta,
                (B0[0] + self.m1 * B1[0] + self.m12 * B2[0])
                    + (B0[1] + self.m1 * B1[1] + self.m12 * B2[1]) * eta
                    + (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * eta.powi(2)
                    + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * eta.powi(3)
                    + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * eta.powi(4)
                    + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * eta.powi(5)
                    + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * eta.powi(6),
            );
        }
        self.i2T0D0.1
    }
    fn i2T0D1(&mut self, eta: f64) -> f64 {
        if eta != self.i2T0D1.0 {
            self.i2T0D1 = (
                eta,
                (B0[1] + self.m1 * B1[1] + self.m12 * B2[1]) * eta
                    + (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * 2.0 * eta.powi(2)
                    + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 3.0 * eta.powi(3)
                    + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 4.0 * eta.powi(4)
                    + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 5.0 * eta.powi(5)
                    + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 6.0 * eta.powi(6),
            )
        }
        self.i2T0D1.1
    }
    fn i2T0D2(&self, eta: f64) -> f64 {
        (B0[2] + self.m1 * B1[2] + self.m12 * B2[2]) * 2.0 * eta.powi(2)
            + (B0[3] + self.m1 * B1[3] + self.m12 * B2[3]) * 6.0 * eta.powi(3)
            + (B0[4] + self.m1 * B1[4] + self.m12 * B2[4]) * 12.0 * eta.powi(4)
            + (B0[5] + self.m1 * B1[5] + self.m12 * B2[5]) * 20.0 * eta.powi(5)
            + (B0[6] + self.m1 * B1[6] + self.m12 * B2[6]) * 30.0 * eta.powi(6)
    }
}
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
