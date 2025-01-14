use crate::algorithms::{brent_zero, shengjin_roots};
use std::f64::consts::SQRT_2;
use std::iter::zip;
const SQRT2ADD1: f64 = SQRT_2 + 1.0;
const SQRT2SUB1: f64 = SQRT_2 - 1.0;
const R_CONST: f64 = 8.314462618;
const AC_COEF: f64 = 0.45724;
const BC_COEF: f64 = 0.07780;
const EPSILON: f64 = f64::EPSILON * 1E8;
pub struct PrMix {
    params: Vec<PrParams>,
}
impl PrMix {
    pub fn new_fluid(crit_t: &[f64], crit_p: &[f64], omega: &[f64]) -> Self {
        let mut params = Vec::with_capacity(omega.len());
        if omega.len() == crit_t.len() && omega.len() == crit_p.len() {
            for i in 0..omega.len() {
                params.push(PrParams {
                    bc: BC_COEF / crit_p[i] * R_CONST * crit_t[i],
                    kappa: 0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i].powi(2),
                    omega1: 7.0 / 3.0 * (omega[i] + 1.0),
                    crit_t: crit_t[i],
                    crit_p: crit_p[i],
                    sqrt_ac: (AC_COEF / crit_p[i]).sqrt() * R_CONST * crit_t[i],
                    sqrt_a: (0.0, 0.0),
                });
            }
        }
        Self { params }
    }
    #[allow(non_snake_case)]
    pub fn fugcoef_exp(&mut self, temp: f64, pres: f64, z: &Vec<f64>, is_liquid: bool) -> Vec<f64> {
        let sum_zi_mul_sqrt_ai = zip(z, &mut self.params)
            .map(|(zi, i)| zi * i.sqrt_a(temp))
            .sum::<f64>();
        let a_mix = zip(z, &mut self.params)
            .map(|(zi, i)| zi * i.sqrt_a(temp) * sum_zi_mul_sqrt_ai)
            .sum::<f64>();
        let b_mix = zip(z, &self.params).map(|(zi, i)| zi * i.b()).sum::<f64>();
        let (A, B) = (
            a_mix * pres / (R_CONST * temp).powi(2),
            b_mix * pres / (R_CONST * temp),
        );
        let (Zv, Zl) = shengjin_roots(
            B - 1.0,
            A - 3.0 * B.powi(2) - 2.0 * B,
            -A * B + B.powi(2) + B.powi(3),
        );
        let Z = if is_liquid && Zl > 0.0 { Zl } else { Zv };
        self.params
            .iter_mut()
            .map(|i| {
                (i.b() / b_mix * (Z - 1.0)
                    - (Z - B).ln()
                    - A / (2.0 * SQRT_2 * B)
                        * (2.0 * i.sqrt_a(temp) * sum_zi_mul_sqrt_ai / a_mix - i.b() / b_mix)
                        * ((Z + SQRT2ADD1 * B) / (Z - SQRT2SUB1 * B)).abs().ln())
                .exp()
            })
            .collect()
    }
    pub fn tpz_flash(&mut self, temp: f64, pres: f64, z: &Vec<f64>) {
        let mut alpha: f64;
        let mut k = self
            .params
            .iter()
            .map(|i| i.p_s(temp) / pres)
            .collect::<Vec<f64>>();
        let mut k0 = k.clone();
        for i in 0..100 {
            alpha = brent_alpha(z, &k);
            let x = zip(z, &k)
                .map(|(zi, ki)| zi / (1.0 + alpha * (ki - 1.0)))
                .collect::<Vec<f64>>();
            let y = zip(z, &k)
                .map(|(zi, ki)| ki * zi / (1.0 + alpha * (ki - 1.0)))
                .collect::<Vec<f64>>();
            let fugcoef_exp_l = self.fugcoef_exp(temp, pres, &x, true);
            let fugcoef_exp_v = self.fugcoef_exp(temp, pres, &y, false);
            k = zip(fugcoef_exp_l, fugcoef_exp_v)
                .map(|(l, v)| l / v)
                .collect::<Vec<f64>>();
            if zip(&k, &k0)
                .map(|(ki, k0i)| (ki - k0i).abs())
                .fold(f64::NEG_INFINITY, |x, y| x.max(y))
                < EPSILON
            {
                println!("i={},x={:?}", i, x);
                println!("i={},y={:?}", i, y);
                println!("i={},k={:?}", i, k);
                break;
            } else {
                k0 = k.clone();
            }
        }
    }
}
struct PrParams {
    bc: f64,
    kappa: f64,
    omega1: f64,
    crit_t: f64,
    crit_p: f64,
    sqrt_ac: f64,
    sqrt_a: (f64, f64),
}
impl PrParams {
    #[inline]
    fn b(&self) -> f64 {
        self.bc
    }
    fn p_s(&self, temp: f64) -> f64 {
        self.crit_p * 10_f64.powf(self.omega1 * (1.0 - self.crit_t / temp))
    }
    fn sqrt_a(&mut self, temp: f64) -> f64 {
        if temp != self.sqrt_a.0 {
            self.sqrt_a = (
                temp,
                self.sqrt_ac * (1.0 + self.kappa * (1.0 - (temp / self.crit_t).sqrt())),
            )
        }
        self.sqrt_a.1
    }
}
/// The mid function
fn brent_alpha(z: &[f64], k: &[f64]) -> f64 {
    let alpha = brent_zero(
        |alpha| {
            -zip(z, k)
                .map(|(zi, ki)| zi * (ki - 1.0) / (1.0 + alpha * (ki - 1.0)))
                .sum::<f64>()
        },
        0.0,
        1.0,
    );
    if alpha.is_nan() {
        if zip(z, k)
            .map(|(zi, ki)| 2.0 * zi * (ki - 1.0) / (ki + 1.0))
            .sum::<f64>()
            .is_sign_negative()
        {
            0.0 // is_sign_positive()
        } else {
            1.0 // is_sign_negative()
        }
    } else {
        alpha
    }
}
/// unit test
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pr_mix() {
        let n2 = (0.025, 126.19, 3.3958e6, 0.03720);
        let c1 = (0.650, 190.56, 4.5992e6, 0.01142);
        let c2 = (0.150, 305.32, 4.8722e6, 0.09950);
        let c3 = (0.150, 369.89, 4.2512e6, 0.15210);
        let c4 = (0.025, 425.13, 3.7960e6, 0.20100);
        let mut example1 = PrMix::new_fluid(
            &[n2.1, c1.1, c2.1, c3.1, c4.1],
            &[n2.2, c1.2, c2.2, c3.2, c4.2],
            &[n2.3, c1.3, c2.3, c3.3, c4.3],
        );
        example1.tpz_flash(250.0, 5.5e6, &vec![n2.0, c1.0, c2.0, c3.0, c4.0]);
    }
}
