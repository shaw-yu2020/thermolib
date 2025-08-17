use crate::f64consts::K;
use std::f64::consts::PI;
const PI2: f64 = PI * PI;
use std::iter::zip;
#[derive(Clone)]
pub struct PolarTerm {
    index_qq: Vec<usize>,
    index_dd: Vec<usize>,
    qq: Option<PqqTerm>,
    dd: Option<PddTerm>,
    dq: Option<PdqTerm>,
}
impl PolarTerm {
    pub fn new(x: &[f64], m: &[f64], sigma: &[f64], epsilon: &[f64], p: &[f64], n: &[i32]) -> Self {
        // qq
        let index_qq: Vec<usize> = n
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n == 4)
            .map(|(i, _)| i)
            .collect();
        let x_qq: Vec<f64> = index_qq.iter().map(|&i| x[i]).collect();
        let m_qq: Vec<f64> = index_qq.iter().map(|&i| m[i]).collect();
        let sigma_qq: Vec<f64> = index_qq.iter().map(|&i| sigma[i]).collect();
        let epsilon_qq: Vec<f64> = index_qq.iter().map(|&i| epsilon[i]).collect();
        let q: Vec<f64> = index_qq.iter().map(|&i| p[i]).collect();
        let qq = if index_qq.is_empty() {
            None
        } else {
            Some(PqqTerm::new(&x_qq, &m_qq, &sigma_qq, &epsilon_qq, &q))
        };
        // dd
        let index_dd: Vec<usize> = n
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n == 2)
            .map(|(i, _)| i)
            .collect();
        let x_dd: Vec<f64> = index_dd.iter().map(|&i| x[i]).collect();
        let m_dd: Vec<f64> = index_dd.iter().map(|&i| m[i]).collect();
        let sigma_dd: Vec<f64> = index_dd.iter().map(|&i| sigma[i]).collect();
        let epsilon_dd: Vec<f64> = index_dd.iter().map(|&i| epsilon[i]).collect();
        let mu: Vec<f64> = index_dd.iter().map(|&i| p[i]).collect();
        let dd = if index_dd.is_empty() {
            None
        } else {
            Some(PddTerm::new(&x_dd, &m_dd, &sigma_dd, &epsilon_dd, &mu))
        };
        // dq
        Self {
            qq,
            dd,
            dq: if index_dd.is_empty() || index_qq.is_empty() {
                None
            } else {
                Some(PdqTerm::new(
                    (&x_dd, &x_qq),
                    (&m_dd, &m_qq),
                    (&sigma_dd, &sigma_qq),
                    (&epsilon_dd, &epsilon_qq),
                    (&mu, &q),
                ))
            },
            index_qq,
            index_dd,
        }
    }
    pub fn new_fracs(&self, x: &[f64]) -> Self {
        let x_qq: Vec<f64> = self.index_qq.iter().map(|&i| x[i]).collect();
        let x_dd: Vec<f64> = self.index_dd.iter().map(|&i| x[i]).collect();
        Self {
            qq: self.qq.as_ref().map(|qq| qq.new_fracs(&x_qq)),
            dd: self.dd.as_ref().map(|dd| dd.new_fracs(&x_dd)),
            dq: self.dq.as_ref().map(|dq| dq.new_fracs((&x_dd, &x_qq))),
            ..self.clone()
        }
    }
    pub fn mu_k<'a>(
        &mut self,
        temp: f64,
        rho_num: f64,
        eta: f64,
        eta_k: &'a [f64],
    ) -> impl Iterator<Item = f64> + use<'a> {
        let coef_eta = rho_num
            * (self
                .qq
                .as_mut()
                .map_or(0.0, |qq| qq.eta1(temp, rho_num, eta))
                + self
                    .dd
                    .as_mut()
                    .map_or(0.0, |dd| dd.eta1(temp, rho_num, eta))
                + self
                    .dq
                    .as_mut()
                    .map_or(0.0, |dq| dq.eta1(temp, rho_num, eta)));
        let mut mu_k: Vec<f64> = eta_k.iter().map(|eta_k| coef_eta * eta_k).collect();
        if let Some(qq) = self.qq.as_mut() {
            zip(qq.mu_qq(temp, rho_num, eta), self.index_qq.iter())
                .map(|(qq, &i)| mu_k[i] += qq)
                .count();
        }
        if let Some(dd) = self.dd.as_mut() {
            zip(dd.mu_dd(temp, rho_num, eta), self.index_dd.iter())
                .map(|(dd, &i)| mu_k[i] += dd)
                .count();
        }
        if let Some(dq) = self.dq.as_mut() {
            zip(dq.mu_dd(temp, rho_num, eta), self.index_dd.iter())
                .map(|(dq, &i)| mu_k[i] += dq)
                .count();
            zip(dq.mu_qq(temp, rho_num, eta), self.index_qq.iter())
                .map(|(dq, &i)| mu_k[i] += dq)
                .count();
        }
        mu_k.into_iter()
    }
    pub fn t0d0(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t0d0(temp, rho_num, eta))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t0d0(temp, rho_num, eta))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t0d0(temp, rho_num, eta))
    }
    pub fn t0d1(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t0d1(temp, rho_num, eta))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t0d1(temp, rho_num, eta))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t0d1(temp, rho_num, eta))
    }
    pub fn t0d2(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t0d2(temp, rho_num, eta))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t0d2(temp, rho_num, eta))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t0d2(temp, rho_num, eta))
    }
    pub fn t0d3(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t0d3(temp, rho_num, eta))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t0d3(temp, rho_num, eta))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t0d3(temp, rho_num, eta))
    }
    pub fn t0d4(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t0d4(temp, rho_num, eta))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t0d4(temp, rho_num, eta))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t0d4(temp, rho_num, eta))
    }
    pub fn t1d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t1d0(temp, rho_num, eta, eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t1d0(temp, rho_num, eta, eta1))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t1d0(temp, rho_num, eta, eta1))
    }
    pub fn t1d1(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t1d1(temp, rho_num, eta, eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t1d1(temp, rho_num, eta, eta1))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t1d1(temp, rho_num, eta, eta1))
    }
    pub fn t1d2(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t1d2(temp, rho_num, eta, eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t1d2(temp, rho_num, eta, eta1))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t1d2(temp, rho_num, eta, eta1))
    }
    pub fn t1d3(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t1d3(temp, rho_num, eta, eta1))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t1d3(temp, rho_num, eta, eta1))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t1d3(temp, rho_num, eta, eta1))
    }
    pub fn t2d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64, eta2: f64) -> f64 {
        self.qq
            .as_mut()
            .map_or(0.0, |qq| qq.t2d0(temp, rho_num, eta, eta1, eta2))
            + self
                .dd
                .as_mut()
                .map_or(0.0, |dd| dd.t2d0(temp, rho_num, eta, eta1, eta2))
            + self
                .dq
                .as_mut()
                .map_or(0.0, |dq| dq.t2d0(temp, rho_num, eta, eta1, eta2))
    }
}
/// macro_rules! fn_polar
macro_rules! fn_polar {
    ($name:ty) => {
        impl $name {
            #[inline]
            fn td_flash(&mut self, temp: f64, rho_num: f64) {
                self.temp = temp;
                self.dens1temp2 = rho_num / (temp * temp);
                self.dens2temp3 = (rho_num * rho_num) / (temp * temp * temp);
            }
            fn eta1(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2eta1(eta) * self.a2t0d0(eta) * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
                    + self.a3eta1(eta) * self.a2t0d0(eta).powi(2))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2)
            }
            fn t0d0(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.td_flash(temp, rho_num);
                self.a2t0d0(eta).powi(2) / (self.a2t0d0(eta) - self.a3t0d0(eta))
            }
            fn t0d1(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t0d1(eta) * self.a2t0d0(eta) * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
                    + self.a3t0d1(eta) * self.a2t0d0(eta).powi(2))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2)
            }
            fn t0d2(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t0d2(eta)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(2) - 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                        + 2.0 * self.a3t0d0(eta).powi(2))
                    + self.a3t0d2(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0 * self.a2t0d1(eta).powi(2) * self.a3t0d0(eta).powi(2)
                    + 2.0 * self.a3t0d1(eta).powi(2) * self.a2t0d0(eta).powi(2)
                    - 4.0
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(3)
            }
            fn t0d3(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t0d3(eta)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(3)
                        - 4.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                        + 5.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                        - 2.0 * self.a3t0d0(eta).powi(3))
                    + self.a3t0d3(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a2t0d2(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a3t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0 * self.a2t0d1(eta).powi(3) * self.a3t0d0(eta).powi(2)
                    + 6.0
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta) + self.a3t0d0(eta))
                    + 6.0 * self.a3t0d1(eta).powi(3) * self.a2t0d0(eta).powi(2)
                    - 6.0
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta) + self.a2t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(4)
            }
            fn t0d4(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t0d4(eta)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(4)
                        - 5.0 * self.a2t0d0(eta).powi(3) * self.a3t0d0(eta)
                        + 9.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta).powi(2)
                        - 7.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(3)
                        + 2.0 * self.a3t0d0(eta).powi(4))
                    + self.a3t0d4(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(3)
                            - 3.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                            + 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                            - self.a3t0d0(eta).powi(3))
                    + 8.0
                        * self.a2t0d3(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 8.0
                        * self.a2t0d3(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 8.0
                        * self.a3t0d3(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 8.0
                        * self.a3t0d3(eta)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 12.0
                        * self.a2t0d2(eta)
                        * self.a3t0d2(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 36.0
                        * self.a2t0d2(eta)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 12.0
                        * self.a2t0d2(eta)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a3t0d0(eta).powi(2))
                    + 24.0
                        * self.a2t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 36.0
                        * self.a3t0d2(eta)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 12.0
                        * self.a3t0d2(eta)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta)
                        * (self.a3t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a2t0d0(eta).powi(2))
                    + 24.0
                        * self.a3t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 24.0 * self.a2t0d1(eta).powi(4) * self.a3t0d0(eta).powi(2)
                    + 24.0 * self.a3t0d1(eta).powi(4) * self.a2t0d0(eta).powi(2)
                    + 24.0
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2)
                            + 4.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 48.0
                        * self.a2t0d1(eta).powi(3)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 48.0
                        * self.a3t0d1(eta).powi(3)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(5)
            }
            fn t1d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t1d0(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
                    + self.a3t1d0(eta, eta1) * self.a2t0d0(eta).powi(2))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2)
            }
            fn t1d1(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t1d1(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(2) - 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                        + 2.0 * self.a3t0d0(eta).powi(2))
                    + self.a3t1d1(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0 * self.a2t1d0(eta, eta1) * self.a2t0d1(eta) * self.a3t0d0(eta).powi(2)
                    - 2.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                    + 2.0 * self.a3t1d0(eta, eta1) * self.a3t0d1(eta) * self.a2t0d0(eta).powi(2)
                    - 2.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(3)
            }
            fn t1d2(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t1d2(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(3)
                        - 4.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                        + 5.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                        - 2.0 * self.a3t0d0(eta).powi(3))
                    + self.a3t1d2(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 4.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 4.0
                        * self.a2t1d1(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 4.0
                        * self.a3t1d1(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 4.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0
                        * self.a2t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 2.0
                        * self.a2t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0
                        * self.a3t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 2.0
                        * self.a3t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                    - 2.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + 2.0 * self.a3t0d0(eta))
                    + 4.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta) + self.a3t0d0(eta))
                    + 6.0
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                    + 2.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 4.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta) + self.a2t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(4)
            }
            fn t1d3(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t1d3(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(4)
                        - 5.0 * self.a2t0d0(eta).powi(3) * self.a3t0d0(eta)
                        + 9.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta).powi(2)
                        - 7.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(3)
                        + 2.0 * self.a3t0d0(eta).powi(4))
                    + self.a3t1d3(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(3)
                            - 3.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                            + 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                            - self.a3t0d0(eta).powi(3))
                    + 6.0
                        * self.a2t1d2(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t1d2(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a2t1d2(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a3t1d2(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 2.0
                        * self.a2t0d3(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 2.0
                        * self.a3t0d3(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 2.0
                        * self.a2t0d3(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 2.0
                        * self.a3t0d3(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d2(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a2t1d1(eta, eta1)
                        * self.a3t0d2(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t1d1(eta, eta1)
                        * self.a3t0d2(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d2(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 18.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 18.0
                        * self.a3t1d1(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a2t1d1(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta)
                        * (self.a3t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a2t0d0(eta).powi(2))
                    + 12.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 12.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    - 18.0
                        * self.a2t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 18.0
                        * self.a3t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 24.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(3)
                        * self.a3t0d0(eta).powi(2)
                    + 24.0
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(3)
                        * self.a2t0d0(eta).powi(2)
                    - 36.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 36.0
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    + 12.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d1(eta)
                        * (self.a2t0d0(eta).powi(2)
                            + 4.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 12.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta)
                        * (self.a2t0d0(eta).powi(2)
                            + 4.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 12.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(3)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 12.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(3)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(5)
            }
            fn t2d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64, eta2: f64) -> f64 {
                self.td_flash(temp, rho_num);
                (self.a2t2d0(eta, eta1, eta2)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(2) - 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                        + 2.0 * self.a3t0d0(eta).powi(2))
                    + self.a3t2d0(eta, eta1, eta2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0 * self.a2t1d0(eta, eta1).powi(2) * self.a3t0d0(eta).powi(2)
                    + 2.0 * self.a3t1d0(eta, eta1).powi(2) * self.a2t0d0(eta).powi(2)
                    - 4.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(3)
            }
            fn a2eta1(&mut self, eta: f64) -> f64 {
                self.vec_a2_coef
                    .iter()
                    .zip(self.vec_epsilon_ij.iter())
                    .zip(self.vec_j2a.iter_mut())
                    .zip(self.vec_j2b.iter_mut())
                    .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                        a2_coef * (j2a.eta1(eta) + epsilon_ij / self.temp * j2b.eta1(eta))
                    })
                    .sum::<f64>()
                    * self.dens1temp2
            }
            fn a2t0d0(&mut self, eta: f64) -> f64 {
                if eta != self.a2t0d0.0 {
                    self.a2t0d0 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef * (j2a.t0d0(eta) + epsilon_ij / self.temp * j2b.t0d0(eta))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t0d0.1
            }
            fn a2t0d1(&mut self, eta: f64) -> f64 {
                if eta != self.a2t0d1.0 {
                    self.a2t0d1 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef * (j2a.t0d1(eta) + epsilon_ij / self.temp * j2b.t0d1(eta))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t0d1.1
            }
            fn a2t0d2(&mut self, eta: f64) -> f64 {
                if eta != self.a2t0d2.0 {
                    self.a2t0d2 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef * (j2a.t0d2(eta) + epsilon_ij / self.temp * j2b.t0d2(eta))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t0d2.1
            }
            fn a2t0d3(&mut self, eta: f64) -> f64 {
                if eta != self.a2t0d3.0 {
                    self.a2t0d3 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef * (j2a.t0d3(eta) + epsilon_ij / self.temp * j2b.t0d3(eta))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t0d3.1
            }
            fn a2t0d4(&mut self, eta: f64) -> f64 {
                if eta != self.a2t0d4.0 {
                    self.a2t0d4 = (
                        eta,
                        /*
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef * (j2a.t0d4(eta) + epsilon_ij / self.temp * j2b.t0d4(eta))
                            })
                            .sum::<f64>()
                            * self.dens1temp2
                         */
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_j2a.iter_mut())
                            .map(|(a2_coef, j2a)| a2_coef * j2a.t0d4(eta))
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t0d4.1
            }
            fn a2t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a2t1d0.0 {
                    self.a2t1d0 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef
                                    * (j2a.t1d0(eta, eta1)
                                        + epsilon_ij / self.temp * j2b.t1d0(eta, eta1))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t1d0.1
            }
            fn a2t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a2t1d1.0 {
                    self.a2t1d1 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef
                                    * (j2a.t1d1(eta, eta1)
                                        + epsilon_ij / self.temp * j2b.t1d1(eta, eta1))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t1d1.1
            }
            fn a2t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a2t1d2.0 {
                    self.a2t1d2 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef
                                    * (j2a.t1d2(eta, eta1)
                                        + epsilon_ij / self.temp * j2b.t1d2(eta, eta1))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t1d2.1
            }
            fn a2t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a2t1d3.0 {
                    self.a2t1d3 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef
                                    * (j2a.t1d3(eta, eta1)
                                        + epsilon_ij / self.temp * j2b.t1d3(eta, eta1))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t1d3.1
            }
            fn a2t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                if eta != self.a2t2d0.0 {
                    self.a2t2d0 = (
                        eta,
                        self.vec_a2_coef
                            .iter()
                            .zip(self.vec_epsilon_ij.iter())
                            .zip(self.vec_j2a.iter_mut())
                            .zip(self.vec_j2b.iter_mut())
                            .map(|(((a2_coef, epsilon_ij), j2a), j2b)| {
                                a2_coef
                                    * (j2a.t2d0(eta, eta1, eta2)
                                        + epsilon_ij / self.temp * j2b.t2d0(eta, eta1, eta2))
                            })
                            .sum::<f64>()
                            * self.dens1temp2,
                    )
                }
                self.a2t2d0.1
            }
            fn a3eta1(&mut self, eta: f64) -> f64 {
                self.vec_a3_coef
                    .iter()
                    .zip(self.vec_j3c.iter_mut())
                    .map(|(a3_coef, j3c)| a3_coef * j3c.eta1(eta))
                    .sum::<f64>()
                    * self.dens2temp3
            }
            fn a3t0d0(&mut self, eta: f64) -> f64 {
                if eta != self.a3t0d0.0 {
                    self.a3t0d0 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t0d0(eta))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t0d0.1
            }
            fn a3t0d1(&mut self, eta: f64) -> f64 {
                if eta != self.a3t0d1.0 {
                    self.a3t0d1 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t0d1(eta))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t0d1.1
            }
            fn a3t0d2(&mut self, eta: f64) -> f64 {
                if eta != self.a3t0d2.0 {
                    self.a3t0d2 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t0d2(eta))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t0d2.1
            }
            fn a3t0d3(&mut self, eta: f64) -> f64 {
                if eta != self.a3t0d3.0 {
                    self.a3t0d3 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t0d3(eta))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t0d3.1
            }
            fn a3t0d4(&mut self, eta: f64) -> f64 {
                if eta != self.a3t0d4.0 {
                    self.a3t0d4 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t0d4(eta))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t0d4.1
            }
            fn a3t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a3t1d0.0 {
                    self.a3t1d0 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t1d0(eta, eta1))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t1d0.1
            }
            fn a3t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a3t1d1.0 {
                    self.a3t1d1 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t1d1(eta, eta1))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t1d1.1
            }
            fn a3t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a3t1d2.0 {
                    self.a3t1d2 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t1d2(eta, eta1))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t1d2.1
            }
            fn a3t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                if eta != self.a3t1d3.0 {
                    self.a3t1d3 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t1d3(eta, eta1))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t1d3.1
            }
            fn a3t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                if eta != self.a3t2d0.0 {
                    self.a3t2d0 = (
                        eta,
                        self.vec_a3_coef
                            .iter()
                            .zip(self.vec_j3c.iter_mut())
                            .map(|(a3_coef, j3c)| a3_coef * j3c.t2d0(eta, eta1, eta2))
                            .sum::<f64>()
                            * self.dens2temp3,
                    )
                }
                self.a3t2d0.1
            }
        }
    };
}
/// J2aTerm
#[derive(Clone)]
struct J2aTerm {
    a: Vec<f64>,
    eta: Vec<f64>,
    t0d0: (f64, f64), // cached variables
    t0d1: (f64, f64), // cached variables
    t0d2: (f64, f64), // cached variables
    t0d3: (f64, f64), // cached variables
}
impl J2aTerm {
    fn new<const N: i8>(m: f64) -> Self {
        let m1 = (m - 1.0) / m; // (m-1)/m
        let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
        Self {
            a: match N {
                4 => vec![
                    QQ_A00 + m1 * QQ_A10 + m12 * QQ_A20,
                    QQ_A01 + m1 * QQ_A11 + m12 * QQ_A21,
                    QQ_A02 + m1 * QQ_A12 + m12 * QQ_A22,
                    QQ_A03 + m1 * QQ_A13 + m12 * QQ_A23,
                    QQ_A04 + m1 * QQ_A14 + m12 * QQ_A24,
                ], // QQ
                2 => vec![
                    DD_A00 + m1 * DD_A10 + m12 * DD_A20,
                    DD_A01 + m1 * DD_A11 + m12 * DD_A21,
                    DD_A02 + m1 * DD_A12 + m12 * DD_A22,
                    DD_A03 + m1 * DD_A13 + m12 * DD_A23,
                    DD_A04 + m1 * DD_A14 + m12 * DD_A24,
                ], // DD
                _ => vec![
                    DQ_A00 + m1 * DQ_A10 + m12 * DQ_A20,
                    DQ_A01 + m1 * DQ_A11 + m12 * DQ_A21,
                    DQ_A02 + m1 * DQ_A12 + m12 * DQ_A22,
                    DQ_A03 + m1 * DQ_A13,
                ], // DQ
            },
            eta: vec![0.0, 0.0, 0.0, 0.0, 0.0],
            // cached variables
            t0d0: (0.0, 0.0),
            t0d1: (0.0, 0.0),
            t0d2: (0.0, 0.0),
            t0d3: (0.0, 0.0),
        }
    }
    #[inline]
    fn eta_flash(&mut self, eta: f64) {
        if eta != self.eta[1] {
            self.eta = vec![1.0, eta, eta.powi(2), eta.powi(3), eta.powi(4)]
        }
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        self.a[1..]
            .iter()
            .zip(self.eta[..4].iter())
            .zip([1.0, 2.0, 3.0, 4.0].iter())
            .map(|((a, eta), coef)| a * eta * coef)
            .sum()
    }
    /// equal to = [rho/T^2 *J2]_t0d0 / {rho/T^2}
    /// equal to = J2t0d0
    fn t0d0(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        if eta != self.t0d0.0 {
            self.t0d0 = (
                eta,
                self.a
                    .iter()
                    .zip(self.eta.iter())
                    .map(|(a, eta)| a * eta)
                    .sum(),
            )
        }
        self.t0d0.1
    }
    /// equal to = [rho/T^2 *J2]_t0d1 / {rho/T^2}
    /// equal to = J2t0d0 + rho * J2t0d1
    fn t0d1(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 1 + n
        if eta != self.t0d1.0 {
            self.t0d1 = (
                eta,
                self.a
                    .iter()
                    .zip(self.eta.iter())
                    .zip([1.0, 2.0, 3.0, 4.0, 5.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum(),
            )
        }
        self.t0d1.1
    }
    /// equal to = [rho/T^2 *J2]_t0d2 / {rho/T^2}
    /// equal to = 2 * rho * J2t0d1 + rho^2 * J2t0d2
    fn t0d2(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 2 * n + ( n - 1 ) * n = ( n + 1 ) * n
        if eta != self.t0d2.0 {
            self.t0d2 = (
                eta,
                self.a[1..]
                    .iter()
                    .zip(self.eta[1..].iter())
                    .zip([2.0, 6.0, 12.0, 20.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum(),
            )
        }
        self.t0d2.1
    }
    /// equal to = [rho/T^2 *J2]_t0d3 / {rho/T^2}
    /// equal to = 3 * rho^2 * J2t0d2 + rho^3 * J2t0d3
    fn t0d3(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 3 * ( n - 1 ) * n + ( n - 2 ) * ( n - 1 ) * n = ( n + 1 ) * ( n - 1 ) * n
        if eta != self.t0d3.0 {
            self.t0d3 = (
                eta,
                self.a[2..]
                    .iter()
                    .zip(self.eta[2..].iter())
                    .zip([6.0, 24.0, 60.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum(),
            )
        }
        self.t0d3.1
    }
    /// equal to = [rho/T^2 *J2]_t0d4 / {rho/T^2}
    /// equal to = 4 * rho^3 * J3t0d3 + rho^4 * J2t0d4
    fn t0d4(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 4 * ( n - 2 ) * ( n - 1 ) * n + ( n - 3 ) * ( n - 2 ) * ( n - 1 ) * n
        // = ( n + 1 ) * ( n - 2 ) * ( n - 1 ) * n
        self.a[3..]
            .iter()
            .zip(self.eta[3..].iter())
            .zip([24.0, 120.0].iter())
            .map(|((a, eta), coef)| a * eta * coef)
            .sum()
    }
    /// equal to = [rho/T^2 *J2]_t1d0 / {rho/T^2}
    /// equal to = -2 * J2t0d0 + T * J2t1d0
    fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -2.0 * self.t0d0(eta)
            + eta1
                * self.a[1..]
                    .iter()
                    .zip(self.eta[..4].iter())
                    .zip([1.0, 2.0, 3.0, 4.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho/T^2 *J2]_t1d1 / {rho/T^2}
    /// equal to = -2 * ( J2t0d0 + rho * J2t0d1 )
    ///            +T * ( J2t1d0 + rho * J2t1d1 )
    fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -2.0 * self.t0d1(eta)
            + eta1
                * self.a[1..]
                    .iter()
                    .zip(self.eta[..4].iter())
                    .zip([2.0, 6.0, 12.0, 20.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho/T^2 *J2]_t1d2 / {rho/T^2}
    /// equal to = -2 * ( 2 * rho * J2t0d1 + rho^2 * J2t0d2 )
    ///            +T * ( 2 * rho * J2t1d1 + rho^2 * J2t1d2 )
    fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -2.0 * self.t0d2(eta)
            + eta1
                * self.a[1..]
                    .iter()
                    .zip(self.eta[..4].iter())
                    .zip([2.0, 12.0, 36.0, 80.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho/T^2 *J2]_t1d3 / {rho/T^2}
    /// equal to = -2 * ( 6 * rho * J2t0d1 + 6 * rho^2 * J2t0d1 + rho^3 * J2t0d3 )
    ///            +T * ( 6 * rho * J2t1d1 + 6 * rho^2 * J2t1d2 + rho^3 * J2t1d3 )
    fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -2.0 * self.t0d3(eta)
            + eta1
                * self.a[2..]
                    .iter()
                    .zip(self.eta[1..4].iter())
                    .zip([12.0, 72.0, 240.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho/T^2 *J2]_t2d0 / {rho/T^2}
    /// equal to = 6 * J2t0d0 - 4 * T * J2t1d0 + T^2 * J2t2d0
    fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        self.eta_flash(eta);
        6.0 * self.t0d0(eta)
            + (eta2 - 4.0 * eta1)
                * self.a[1..]
                    .iter()
                    .zip(self.eta[..4].iter())
                    .zip([1.0, 2.0, 3.0, 4.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum::<f64>()
            + eta1.powi(2)
                * self.a[2..]
                    .iter()
                    .zip(self.eta[..3].iter())
                    .zip([2.0, 6.0, 12.0].iter())
                    .map(|((a, eta), coef)| a * eta * coef)
                    .sum::<f64>()
    }
}
const QQ_A00: f64 = 1.2378308;
const QQ_A01: f64 = 2.4355031;
const QQ_A02: f64 = 1.6330905;
const QQ_A03: f64 = -1.6118152;
const QQ_A04: f64 = 6.9771185;
const QQ_A10: f64 = 1.2854109;
const QQ_A11: f64 = -11.465615;
const QQ_A12: f64 = 22.086893;
const QQ_A13: f64 = 7.4691383;
const QQ_A14: f64 = -17.197772;
const QQ_A20: f64 = 1.7942954;
const QQ_A21: f64 = 0.7695103;
const QQ_A22: f64 = 7.2647923;
const QQ_A23: f64 = 94.486699;
const QQ_A24: f64 = -77.148458;
const DD_A00: f64 = 0.3043504;
const DD_A01: f64 = -0.1358588;
const DD_A02: f64 = 1.4493329;
const DD_A03: f64 = 0.3556977;
const DD_A04: f64 = -2.0653308;
const DD_A10: f64 = 0.9534641;
const DD_A11: f64 = -1.8396383;
const DD_A12: f64 = 2.0131180;
const DD_A13: f64 = -7.3724958;
const DD_A14: f64 = 8.2374135;
const DD_A20: f64 = -1.1610080;
const DD_A21: f64 = 4.5258607;
const DD_A22: f64 = 0.9751222;
const DD_A23: f64 = -12.281038;
const DD_A24: f64 = 5.9397575;
const DQ_A00: f64 = 0.6970950;
const DQ_A01: f64 = -0.6335541;
const DQ_A02: f64 = 2.9455090;
const DQ_A03: f64 = -1.4670273;
const DQ_A10: f64 = -0.6734593;
const DQ_A11: f64 = -1.4258991;
const DQ_A12: f64 = 4.1944139;
const DQ_A13: f64 = 1.0266216;
const DQ_A20: f64 = 0.6703408;
const DQ_A21: f64 = -4.3384718;
const DQ_A22: f64 = 7.2341684;
/// J2bTerm
#[derive(Clone)]
struct J2bTerm {
    b0: f64,
    b1: f64,
    b2: f64,
    t0d0: (f64, f64), // cached variables
    t0d1: (f64, f64), // cached variables
    t0d2: (f64, f64), // cached variables
    t0d3: (f64, f64), // cached variables
}
impl J2bTerm {
    fn new<const N: i8>(m: f64) -> Self {
        let m1 = (m - 1.0) / m; // (m-1)/m
        let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
        Self {
            b0: match N {
                4 => QQ_B00 + m1 * QQ_B10 + m12 * QQ_B20, // QQ
                2 => DD_B00 + m1 * DD_B10 + m12 * DD_B20, // DD
                _ => DQ_B00 + m1 * DQ_B10 + m12 * DQ_B20, // DQ
            },
            b1: match N {
                4 => QQ_B01 + m1 * QQ_B11 + m12 * QQ_B21, // QQ
                2 => DD_B01 + m1 * DD_B11 + m12 * DD_B21, // DD
                _ => DQ_B01 + m1 * DQ_B11 + m12 * DQ_B21, // DQ
            },
            b2: match N {
                4 => QQ_B02 + m1 * QQ_B12 + m12 * QQ_B22, // QQ
                2 => DD_B02 + m1 * DD_B12 + m12 * DD_B22, // DD
                _ => DQ_B02 + m1 * DQ_B12,                // DQ
            },
            // cached variables
            t0d0: (0.0, 0.0),
            t0d1: (0.0, 0.0),
            t0d2: (0.0, 0.0),
            t0d3: (0.0, 0.0),
        }
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        self.b1 + self.b2 * eta * 2.0
    }
    /// equal to = [rho/T^3 *J2]_t0d0 / {rho/T^3}
    /// equal to = J2t0d0
    fn t0d0(&mut self, eta: f64) -> f64 {
        if eta != self.t0d0.0 {
            self.t0d0 = (eta, self.b0 + self.b1 * eta + self.b2 * eta.powi(2))
        }
        self.t0d0.1
    }
    /// equal to = [rho/T^3 *J2]_t0d1 / {rho/T^3}
    /// equal to = J2t0d0 + rho * J2t0d1
    fn t0d1(&mut self, eta: f64) -> f64 {
        // 1 + n
        if eta != self.t0d1.0 {
            self.t0d1 = (
                eta,
                self.b0 + self.b1 * 2.0 * eta + self.b2 * 3.0 * eta.powi(2),
            )
        }
        self.t0d1.1
    }
    /// equal to = [rho/T^3 *J2]_t0d2 / {rho/T^3}
    /// equal to = 2 * rho * J2t0d1 + rho^2 * J2t0d2
    fn t0d2(&mut self, eta: f64) -> f64 {
        // 2 * n + ( n - 1 ) * n = ( n + 1 ) * n
        if eta != self.t0d2.0 {
            self.t0d2 = (eta, self.b1 * 2.0 * eta + self.b2 * 6.0 * eta.powi(2))
        }
        self.t0d2.1
    }
    /// equal to = [rho/T^3 *J2]_t0d3 / {rho/T^3}
    /// equal to = 3 * rho^2 * J2t0d2 + rho^3 * J2t0d3
    fn t0d3(&mut self, eta: f64) -> f64 {
        // 3 * ( n - 1 ) * n + ( n - 2 ) * ( n - 1 ) * n = ( n + 1 ) * ( n - 1 ) * n
        if eta != self.t0d3.0 {
            self.t0d3 = (eta, self.b2 * 6.0 * eta.powi(2))
        }
        self.t0d3.1
    }
    /*
    /// equal to = [rho/T^3 *J2]_t0d4 / {rho/T^3}
    /// equal to = 4 * rho^3 * J3t0d3 + rho^4 * J2t0d4
    #[inline]
    fn t0d4(&mut self, _eta: f64) -> f64 {
        // 4 * ( n - 2 ) * ( n - 1 ) * n + ( n - 3 ) * ( n - 2 ) * ( n - 1 ) * n
        0.0
    }
     */
    /// equal to = [rho/T^3 *J2]_t1d0 / {rho/T^3}
    /// equal to = -3 * J2t0d0 + T * J2t1d0
    fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        -3.0 * self.t0d0(eta) + eta1 * (self.b1 + self.b2 * 2.0 * eta)
    }
    /// equal to = [rho/T^3 *J2]_t1d1 / {rho/T^3}
    /// equal to = -3 * ( J2t0d0 + rho * J2t0d1 )
    ///            +T * ( J2t1d0 + rho * J2t1d1 )
    fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        -3.0 * self.t0d1(eta) + eta1 * (self.b1 * 2.0 + self.b2 * 6.0 * eta)
    }
    /// equal to = [rho/T^3 *J2]_t1d2 / {rho/T^3}
    /// equal to = -3 * ( 2 * rho * J2t0d1 + rho^2 * J2t0d2 )
    ///            +T * ( 2 * rho * J2t1d1 + rho^2 * J2t1d2 )
    fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        -3.0 * self.t0d2(eta) + eta1 * (self.b1 * 2.0 + self.b2 * 12.0 * eta)
    }
    /// equal to = [rho/T^3 *J2]_t1d3 / {rho/T^3}
    /// equal to = -3 * ( 6 * rho * J2t0d1 + 6 * rho^2 * J2t0d1 + rho^3 * J2t0d3 )
    ///            +T * ( 6 * rho * J2t1d1 + 6 * rho^2 * J2t1d2 + rho^3 * J2t1d3 )
    fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        -3.0 * self.t0d3(eta) + eta1 * self.b2 * 12.0 * eta
    }
    /// equal to = [rho/T^3 *J2]_t2d0 / {rho/T^3}
    /// equal to = 12 * J2t0d0 - 6 * T * J2t1d0 + T^2 * J2t2d0
    fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        12.0 * self.t0d0(eta)
            + (eta2 - 6.0 * eta1) * (self.b1 + self.b2 * 2.0 * eta)
            + eta1.powi(2) * self.b2 * 2.0
    }
}
const QQ_B00: f64 = 0.4542718;
const QQ_B01: f64 = -4.5016264;
const QQ_B02: f64 = 3.5858868;
const QQ_B10: f64 = -0.8137340;
const QQ_B11: f64 = 10.064030;
const QQ_B12: f64 = -10.876631;
const QQ_B20: f64 = 6.8682675;
const QQ_B21: f64 = -5.1732238;
const QQ_B22: f64 = -17.240207;
const DD_B00: f64 = 0.2187939;
const DD_B01: f64 = -1.1896431;
const DD_B02: f64 = 1.1626889;
const DD_B10: f64 = -0.5873164;
const DD_B11: f64 = 1.2489132;
const DD_B12: f64 = -0.5085280;
const DD_B20: f64 = 3.4869576;
const DD_B21: f64 = -14.915974;
const DD_B22: f64 = 15.372022;
const DQ_B00: f64 = -0.4840383;
const DQ_B01: f64 = 1.9704055;
const DQ_B02: f64 = -2.1185727;
const DQ_B10: f64 = 0.6765101;
const DQ_B11: f64 = -3.0138675;
const DQ_B12: f64 = 0.4674266;
const DQ_B20: f64 = -1.1675601;
const DQ_B21: f64 = 2.1348843;
/// J3cTerm
#[derive(Clone)]
struct J3cTerm {
    c: Vec<f64>,
    eta: Vec<f64>,
    t0d0: (f64, f64), // cached variables
    t0d1: (f64, f64), // cached variables
    t0d2: (f64, f64), // cached variables
    t0d3: (f64, f64), // cached variables
}
impl J3cTerm {
    fn new<const N: i8>(m: f64) -> Self {
        let m1 = (m - 1.0) / m; // (m-1)/m
        let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
        Self {
            c: match N {
                4 => vec![
                    QQ_C00 + m1 * QQ_C10 + m12 * QQ_C20,
                    QQ_C01 + m1 * QQ_C11 + m12 * QQ_C21,
                    QQ_C02 + m1 * QQ_C12 + m12 * QQ_C22,
                    QQ_C03 + m1 * QQ_C13,
                ], // QQ
                2 => vec![
                    DD_C00 + m1 * DD_C10 + m12 * DD_C20,
                    DD_C01 + m1 * DD_C11 + m12 * DD_C21,
                    DD_C02 + m1 * DD_C12 + m12 * DD_C22,
                    DD_C03 + m1 * DD_C13 + m12 * DD_C23,
                ], // DD
                _ => vec![
                    DQ_C00 + m1 * DQ_C10,
                    DQ_C01 + m1 * DQ_C11,
                    DQ_C02 + m1 * DQ_C12,
                ], // DQ
            },
            eta: vec![0.0, 0.0, 0.0, 0.0],
            // cached variables
            t0d0: (0.0, 0.0),
            t0d1: (0.0, 0.0),
            t0d2: (0.0, 0.0),
            t0d3: (0.0, 0.0),
        }
    }
    #[inline]
    fn eta_flash(&mut self, eta: f64) {
        if eta != self.eta[1] {
            self.eta = vec![1.0, eta, eta.powi(2), eta.powi(3)]
        }
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        self.c[1..]
            .iter()
            .zip(self.eta[..3].iter())
            .zip([1.0, 2.0, 3.0].iter())
            .map(|((c, eta), coef)| c * eta * coef)
            .sum()
    }
    /// equal to = [rho^2/T^3 *J3]_t0d0 / {rho^2/T^3}
    /// equal to = J3t0d0
    fn t0d0(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        if eta != self.t0d0.0 {
            self.t0d0 = (
                eta,
                self.c
                    .iter()
                    .zip(self.eta.iter())
                    .map(|(c, eta)| c * eta)
                    .sum(),
            )
        }
        self.t0d0.1
    }
    /// equal to = [rho^2/T^3 *J3]_t0d1 / {rho^2/T^3}
    /// equal to = 2 * J3t0d0 + rho * J3t0d1
    fn t0d1(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 2 + n
        if eta != self.t0d1.0 {
            self.t0d1 = (
                eta,
                self.c
                    .iter()
                    .zip(self.eta.iter())
                    .zip([2.0, 3.0, 4.0, 5.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum(),
            )
        }
        self.t0d1.1
    }
    /// equal to = [rho^2/T^3 *J3]_t0d2 / {rho^2/T^3}
    /// equal to = 2 * J3t0d0 + 4 * rho * J3t0d1 + rho^2 * J3t0d2
    fn t0d2(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 2 + 4 * n + ( n - 1 ) * n
        // = 2 + 3 * n + n^2
        if eta != self.t0d2.0 {
            self.t0d2 = (
                eta,
                self.c
                    .iter()
                    .zip(self.eta.iter())
                    .zip([2.0, 6.0, 12.0, 20.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum(),
            )
        }
        self.t0d2.1
    }
    /// equal to = [rho^2/T^3 *J3]_t0d3 / {rho^2/T^3}
    /// equal to = rho * ( 6 * J3t0d1 + 6 * rho * J3t0d2 + rho^2 * J3t0d3)
    fn t0d3(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 6 * n + 6 * ( n - 1 ) * n + ( n - 2 ) * ( n - 1 ) * n
        // = ( 2 + 3 * n + n^2 ) * n
        if eta != self.t0d3.0 {
            self.t0d3 = (
                eta,
                self.c[1..]
                    .iter()
                    .zip(self.eta[1..].iter())
                    .zip([6.0, 24.0, 60.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum(),
            )
        }
        self.t0d3.1
    }
    /// equal to = [rho^2/T^3 *J3]_t0d4 / {rho^2/T^3}
    /// equal to = rho^2 * ( 12 * J3t0d2 + 8 * rho * J3t0d3 + rho^2 * J3t0d4)
    fn t0d4(&mut self, eta: f64) -> f64 {
        self.eta_flash(eta);
        // 12 * ( n - 1 ) * n
        // + 8 * ( n - 2 ) * ( n - 1 ) * n
        // + ( n - 3 ) * ( n - 2 ) * ( n - 1 ) * n
        // = ( 2 + 3 * n + n^2 ) * ( n - 1 ) * n
        self.c[2..]
            .iter()
            .zip(self.eta[2..].iter())
            .zip([24.0, 120.0].iter())
            .map(|((c, eta), coef)| c * eta * coef)
            .sum()
    }
    /// equal to = [rho^2/T^3 *J3]_t1d0 / {rho^2/T^3}
    /// equal to = -3 * J3t0d0 + T * J3t1d0
    fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -3.0 * self.t0d0(eta)
            + eta1
                * self.c[1..]
                    .iter()
                    .zip(self.eta[..3].iter())
                    .zip([1.0, 2.0, 3.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho^2/T^3 *J3]_t1d1 / {rho^2/T^3}
    /// equal to = -3 * ( 2 * J3t0d0 + rho * J3t0d1 )
    ///            +T * ( 2 * J3t1d0 + rho * J3t1d1 )
    fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -3.0 * self.t0d1(eta)
            + eta1
                * self.c[1..]
                    .iter()
                    .zip(self.eta[..3].iter())
                    .zip([3.0, 8.0, 15.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho^2/T^3 *J3]_t1d2 / {rho^2/T^3}
    /// equal to = -3 * ( 2 * J3t0d0 + 4 * rho * J3t0d1 + rho^2 * J3t0d2 )
    ///            +T * ( 2 * J3t1d0 + 4 * rho * J3t1d1 + rho^2 * J3t1d2 )
    fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -3.0 * self.t0d2(eta)
            + eta1
                * self.c[1..]
                    .iter()
                    .zip(self.eta[..3].iter())
                    .zip([6.0, 24.0, 60.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho^2/T^3 *J3]_t1d3 / {rho^2/T^3}
    /// equal to = -3 * ( 6 * rho * J3t0d1 + 6 * rho^2 * J3t0d2 + rho^3 * J3t0d3 )
    ///            +T * ( 6 * rho * J3t1d1 + 6 * rho^2 * J3t1d2 + rho^3 * J3t1d3 )
    fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        self.eta_flash(eta);
        -3.0 * self.t0d3(eta)
            + eta1
                * self.c[1..]
                    .iter()
                    .zip(self.eta[..3].iter())
                    .zip([6.0, 48.0, 180.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum::<f64>()
    }
    /// equal to = [rho^2/T^3 *J3]_t2d0 / {rho^2/T^3}
    /// equal to = 12 * J3t0d0 - 6 * T * J3t1d0 + T^2 * J3t2d0
    fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        self.eta_flash(eta);
        12.0 * self.t0d0(eta)
            + (eta2 - 6.0 * eta1)
                * self.c[1..]
                    .iter()
                    .zip(self.eta[..3].iter())
                    .zip([1.0, 2.0, 3.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum::<f64>()
            + eta1.powi(2)
                * self.c[2..]
                    .iter()
                    .zip(self.eta[..2].iter())
                    .zip([2.0, 6.0].iter())
                    .map(|((c, eta), coef)| c * eta * coef)
                    .sum::<f64>()
    }
}
const QQ_C00: f64 = -0.5000437;
const QQ_C01: f64 = 6.5318692;
const QQ_C02: f64 = -16.014780;
const QQ_C03: f64 = 14.425970;
const QQ_C10: f64 = 2.0002094;
const QQ_C11: f64 = -6.7838658;
const QQ_C12: f64 = 20.383246;
const QQ_C13: f64 = -10.895984;
const QQ_C20: f64 = 3.1358271;
const QQ_C21: f64 = 7.2475888;
const QQ_C22: f64 = 3.0759478;
const DD_C00: f64 = -0.0646774;
const DD_C01: f64 = 0.1975882;
const DD_C02: f64 = -0.8087562;
const DD_C03: f64 = 0.6902849;
const DD_C10: f64 = -0.9520876;
const DD_C11: f64 = 2.9924258;
const DD_C12: f64 = -2.3802636;
const DD_C13: f64 = -0.2701261;
const DD_C20: f64 = -0.6260979;
const DD_C21: f64 = 1.2924686;
const DD_C22: f64 = 1.6542783;
const DD_C23: f64 = -3.4396744;
const DQ_C00: f64 = 7.846431;
const DQ_C01: f64 = 33.42700;
const DQ_C02: f64 = 4.689111;
const DQ_C10: f64 = -20.72202;
const DQ_C11: f64 = -58.63904;
const DQ_C12: f64 = -1.764887;
/// PqqTerm
#[derive(Clone)]
struct PqqTerm {
    temp: f64,
    dens1temp2: f64, // dens^1 / temp^2
    dens2temp3: f64, // dens^2 / temp^3
    // a2term
    vec_epsilon_ij: Vec<f64>,
    vec_a2_coef: Vec<f64>,
    vec_j2a: Vec<J2aTerm>,
    vec_j2b: Vec<J2bTerm>,
    // a3term
    vec_a3_coef: Vec<f64>,
    vec_j3c: Vec<J3cTerm>,
    // calc mu_k
    x: Vec<f64>,
    m: Vec<f64>,
    epsilon: Vec<f64>, // for mat_epsilon_ij
    mat_a2_coef: Vec<Vec<f64>>,
    mat_a3_coef: Vec<Vec<Vec<f64>>>,
    a2t0d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d4: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t2d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d4: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t2d0: (f64, f64), // cached variables -> macro_rules! fn_polar
}
fn_polar!(PqqTerm);
impl PqqTerm {
    fn new(x: &[f64], m: &[f64], sigma: &[f64], epsilon: &[f64], q: &[f64]) -> Self {
        let q2_plus: Vec<f64> = q
            .iter()
            .zip(m)
            .map(|(q, m)| q.powi(2) * 1E-19 / m / K)
            .collect();
        let m: Vec<f64> = m.iter().map(|m| m.min(2.0)).collect(); // HHH
        let mat_a2_coef: Vec<Vec<f64>> = q2_plus
            .iter()
            .zip(sigma)
            .map(|(q2_i, sigma_i)| {
                q2_plus
                    .iter()
                    .zip(sigma)
                    .map(|(q2_j, sigma_j)| -72.0 * PI * q2_i * q2_j / (sigma_i + sigma_j).powi(7))
                    .collect()
            })
            .collect();
        let mat_a3_coef: Vec<Vec<Vec<f64>>> = q2_plus
            .iter()
            .zip(sigma)
            .map(|(q2_i, sigma_i)| {
                q2_plus
                    .iter()
                    .zip(sigma)
                    .map(|(q2_j, sigma_j)| {
                        q2_plus
                            .iter()
                            .zip(sigma)
                            .map(|(q2_k, sigma_k)| {
                                288.0 * PI2 * q2_i * q2_j * q2_k
                                    / ((sigma_i + sigma_j)
                                        * (sigma_i + sigma_k)
                                        * (sigma_j + sigma_k))
                                        .powi(3)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let mut vec_epsilon_ij: Vec<f64> = Vec::new();
        let mut vec_a2_coef: Vec<f64> = Vec::new();
        let mut vec_j2a: Vec<J2aTerm> = Vec::new();
        let mut vec_j2b: Vec<J2bTerm> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        vec_epsilon_ij.push((epsilon[i] * epsilon[j]).sqrt());
                        vec_a2_coef.push(2.0 * x[i] * x[j] * mat_a2_coef[i][j]);
                        let m_ij = (m[i] * m[j]).sqrt();
                        vec_j2a.push(J2aTerm::new::<4>(m_ij));
                        vec_j2b.push(J2bTerm::new::<4>(m_ij));
                    })
                    .count();
                vec_epsilon_ij.push(epsilon[i]);
                vec_a2_coef.push(x[i] * x[i] * mat_a2_coef[i][i]);
                vec_j2a.push(J2aTerm::new::<4>(m[i]));
                vec_j2b.push(J2bTerm::new::<4>(m[i]));
            })
            .count();
        let mut vec_a3_coef: Vec<f64> = Vec::new();
        let mut vec_j3c: Vec<J3cTerm> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                vec_a3_coef.push(6.0 * x[i] * x[j] * x[k] * mat_a3_coef[i][j][k]);
                                vec_j3c.push(J3cTerm::new::<4>((m[i] * m[j] * m[k]).cbrt()));
                            })
                            .count();
                        vec_a3_coef.push(3.0 * x[i] * x[i] * x[j] * mat_a3_coef[i][i][j]);
                        vec_j3c.push(J3cTerm::new::<4>((m[i] * m[i] * m[j]).cbrt()));
                        vec_a3_coef.push(3.0 * x[i] * x[j] * x[j] * mat_a3_coef[i][j][j]);
                        vec_j3c.push(J3cTerm::new::<4>((m[i] * m[j] * m[j]).cbrt()));
                    })
                    .count();
                vec_a3_coef.push(x[i].powi(3) * mat_a3_coef[i][i][i]);
                vec_j3c.push(J3cTerm::new::<4>(m[i]));
            })
            .count();
        Self {
            temp: 0.0,
            dens1temp2: 0.0,
            dens2temp3: 0.0,
            // a2term
            vec_epsilon_ij,
            vec_a2_coef,
            vec_j2a,
            vec_j2b,
            // a3term
            vec_a3_coef,
            vec_j3c,
            // calc mu_k
            x: x.to_owned(),
            m: m.to_owned(),
            epsilon: epsilon.to_owned(),
            mat_a2_coef,
            mat_a3_coef,
            a2t0d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d4: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t2d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d4: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t2d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
        }
    }
    fn new_fracs(&self, x: &[f64]) -> Self {
        let mut vec_a2_coef: Vec<f64> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        vec_a2_coef.push(2.0 * x[i] * x[j] * self.mat_a2_coef[i][j]);
                    })
                    .count();
                vec_a2_coef.push(x[i] * x[i] * self.mat_a2_coef[i][i]);
            })
            .count();
        let mut vec_a3_coef: Vec<f64> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                vec_a3_coef
                                    .push(6.0 * x[i] * x[j] * x[k] * self.mat_a3_coef[i][j][k]);
                            })
                            .count();
                        vec_a3_coef.push(3.0 * x[i] * x[i] * x[j] * self.mat_a3_coef[i][i][j]);
                        vec_a3_coef.push(3.0 * x[i] * x[j] * x[j] * self.mat_a3_coef[i][j][j]);
                    })
                    .count();
                vec_a3_coef.push(x[i].powi(3) * self.mat_a3_coef[i][i][i]);
            })
            .count();
        Self {
            x: x.to_owned(),
            vec_a2_coef,
            vec_a3_coef,
            ..self.clone()
        }
    }
    fn mu_qq(&mut self, temp: f64, rho_num: f64, eta: f64) -> Vec<f64> {
        self.td_flash(temp, rho_num);
        let coef_a2 = self.a2t0d0(eta) * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
            / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a2: Vec<f64> = self
            .mat_a2_coef
            .iter()
            .enumerate()
            .map(|(i, a2_coef_i)| {
                2.0 * a2_coef_i
                    .iter()
                    .enumerate()
                    .map(|(j, a2_coef_ij)| {
                        self.x[j]
                            * a2_coef_ij
                            * (J2aTerm::new::<4>((self.m[i] * self.m[j]).sqrt()).t0d0(eta)
                                + (self.epsilon[i] * self.epsilon[j]).sqrt() / self.temp
                                    * J2bTerm::new::<4>((self.m[i] * self.m[j]).sqrt()).t0d0(eta))
                    })
                    .sum::<f64>()
                    * self.dens1temp2
            })
            .collect();
        let coef_a3 = self.a2t0d0(eta).powi(2) / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a3: Vec<f64> = self
            .mat_a3_coef
            .iter()
            .enumerate()
            .map(|(i, a3_coef_i)| {
                3.0 * a3_coef_i
                    .iter()
                    .enumerate()
                    .map(|(j, a3_coef_ij)| {
                        self.x[j]
                            * a3_coef_ij
                                .iter()
                                .enumerate()
                                .map(|(k, a3_coef_ijk)| {
                                    self.x[k]
                                        * a3_coef_ijk
                                        * J3cTerm::new::<4>(
                                            (self.m[i] * self.m[j] * self.m[k]).cbrt(),
                                        )
                                        .t0d0(eta)
                                })
                                .sum::<f64>()
                    })
                    .sum::<f64>()
                    * self.dens2temp3
            })
            .collect();
        zip(a2, a3)
            .map(move |(a2, a3)| coef_a2 * a2 + coef_a3 * a3)
            .collect()
    }
}
/// PddTerm
#[derive(Clone)]
struct PddTerm {
    temp: f64,
    dens1temp2: f64, // dens^1 / temp^2
    dens2temp3: f64, // dens^2 / temp^3
    // a2term
    vec_epsilon_ij: Vec<f64>,
    vec_a2_coef: Vec<f64>,
    vec_j2a: Vec<J2aTerm>,
    vec_j2b: Vec<J2bTerm>,
    // a3term
    vec_a3_coef: Vec<f64>,
    vec_j3c: Vec<J3cTerm>,
    // calc mu_k
    x: Vec<f64>,
    m: Vec<f64>,
    epsilon: Vec<f64>, // for mat_epsilon_ij
    mat_a2_coef: Vec<Vec<f64>>,
    mat_a3_coef: Vec<Vec<Vec<f64>>>,
    a2t0d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d4: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t2d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d4: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t2d0: (f64, f64), // cached variables -> macro_rules! fn_polar
}
fn_polar!(PddTerm);
impl PddTerm {
    fn new(x: &[f64], m: &[f64], sigma: &[f64], epsilon: &[f64], mu: &[f64]) -> Self {
        let mu2_plus: Vec<f64> = mu
            .iter()
            .zip(m)
            .map(|(mu, m)| mu.powi(2) * 1E-19 / m / K)
            .collect();
        let m: Vec<f64> = m.iter().map(|m| m.min(2.0)).collect();
        let mat_a2_coef: Vec<Vec<f64>> = mu2_plus
            .iter()
            .zip(sigma)
            .map(|(mu2_i, sigma_i)| {
                mu2_plus
                    .iter()
                    .zip(sigma)
                    .map(|(mu2_j, sigma_j)| -8.0 * PI * mu2_i * mu2_j / (sigma_i + sigma_j).powi(3))
                    .collect()
            })
            .collect();
        let mat_a3_coef: Vec<Vec<Vec<f64>>> = mu2_plus
            .iter()
            .zip(sigma)
            .map(|(mu2_i, sigma_i)| {
                mu2_plus
                    .iter()
                    .zip(sigma)
                    .map(|(mu2_j, sigma_j)| {
                        mu2_plus
                            .iter()
                            .zip(sigma)
                            .map(|(mu2_k, sigma_k)| {
                                -32.0 / 3.0 * PI2 * mu2_i * mu2_j * mu2_k
                                    / ((sigma_i + sigma_j)
                                        * (sigma_i + sigma_k)
                                        * (sigma_j + sigma_k))
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let mut vec_epsilon_ij: Vec<f64> = Vec::new();
        let mut vec_a2_coef: Vec<f64> = Vec::new();
        let mut vec_j2a: Vec<J2aTerm> = Vec::new();
        let mut vec_j2b: Vec<J2bTerm> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        vec_epsilon_ij.push((epsilon[i] * epsilon[j]).sqrt());
                        vec_a2_coef.push(2.0 * x[i] * x[j] * mat_a2_coef[i][j]);
                        let m_ij = (m[i] * m[j]).sqrt();
                        vec_j2a.push(J2aTerm::new::<2>(m_ij));
                        vec_j2b.push(J2bTerm::new::<2>(m_ij));
                    })
                    .count();
                vec_epsilon_ij.push(epsilon[i]);
                vec_a2_coef.push(x[i] * x[i] * mat_a2_coef[i][i]);
                vec_j2a.push(J2aTerm::new::<2>(m[i]));
                vec_j2b.push(J2bTerm::new::<2>(m[i]));
            })
            .count();
        let mut vec_a3_coef: Vec<f64> = Vec::new();
        let mut vec_j3c: Vec<J3cTerm> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                vec_a3_coef.push(6.0 * x[i] * x[j] * x[k] * mat_a3_coef[i][j][k]);
                                vec_j3c.push(J3cTerm::new::<2>((m[i] * m[j] * m[k]).cbrt()));
                            })
                            .count();
                        vec_a3_coef.push(3.0 * x[i] * x[i] * x[j] * mat_a3_coef[i][i][j]);
                        vec_j3c.push(J3cTerm::new::<2>((m[i] * m[i] * m[j]).cbrt()));
                        vec_a3_coef.push(3.0 * x[i] * x[j] * x[j] * mat_a3_coef[i][j][j]);
                        vec_j3c.push(J3cTerm::new::<2>((m[i] * m[j] * m[j]).cbrt()));
                    })
                    .count();
                vec_a3_coef.push(x[i].powi(3) * mat_a3_coef[i][i][i]);
                vec_j3c.push(J3cTerm::new::<2>(m[i]));
            })
            .count();
        Self {
            temp: 0.0,
            dens1temp2: 0.0,
            dens2temp3: 0.0,
            // a2term
            vec_epsilon_ij,
            vec_a2_coef,
            vec_j2a,
            vec_j2b,
            // a3term
            vec_a3_coef,
            vec_j3c,
            // calc mu_k
            x: x.to_owned(),
            m: m.to_owned(),
            epsilon: epsilon.to_owned(),
            mat_a2_coef,
            mat_a3_coef,
            a2t0d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d4: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t2d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d4: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t2d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
        }
    }
    fn new_fracs(&self, x: &[f64]) -> Self {
        let mut vec_a2_coef: Vec<f64> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        vec_a2_coef.push(2.0 * x[i] * x[j] * self.mat_a2_coef[i][j]);
                    })
                    .count();
                vec_a2_coef.push(x[i] * x[i] * self.mat_a2_coef[i][i]);
            })
            .count();
        let mut vec_a3_coef: Vec<f64> = Vec::new();
        let _ = (0..x.len())
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                vec_a3_coef
                                    .push(6.0 * x[i] * x[j] * x[k] * self.mat_a3_coef[i][j][k]);
                            })
                            .count();
                        vec_a3_coef.push(3.0 * x[i] * x[i] * x[j] * self.mat_a3_coef[i][i][j]);
                        vec_a3_coef.push(3.0 * x[i] * x[j] * x[j] * self.mat_a3_coef[i][j][j]);
                    })
                    .count();
                vec_a3_coef.push(x[i].powi(3) * self.mat_a3_coef[i][i][i]);
            })
            .count();
        Self {
            x: x.to_owned(),
            vec_a2_coef,
            vec_a3_coef,
            ..self.clone()
        }
    }
    fn mu_dd(&mut self, temp: f64, rho_num: f64, eta: f64) -> Vec<f64> {
        self.td_flash(temp, rho_num);
        let coef_a2 = self.a2t0d0(eta) * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
            / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a2: Vec<f64> = self
            .mat_a2_coef
            .iter()
            .enumerate()
            .map(|(i, a2_coef_i)| {
                2.0 * a2_coef_i
                    .iter()
                    .enumerate()
                    .map(|(j, a2_coef_ij)| {
                        self.x[j]
                            * a2_coef_ij
                            * (J2aTerm::new::<2>((self.m[i] * self.m[j]).sqrt()).t0d0(eta)
                                + (self.epsilon[i] * self.epsilon[j]).sqrt() / self.temp
                                    * J2bTerm::new::<2>((self.m[i] * self.m[j]).sqrt()).t0d0(eta))
                    })
                    .sum::<f64>()
                    * self.dens1temp2
            })
            .collect();
        let coef_a3 = self.a2t0d0(eta).powi(2) / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a3: Vec<f64> = self
            .mat_a3_coef
            .iter()
            .enumerate()
            .map(|(i, a3_coef_i)| {
                3.0 * a3_coef_i
                    .iter()
                    .enumerate()
                    .map(|(j, a3_coef_ij)| {
                        self.x[j]
                            * a3_coef_ij
                                .iter()
                                .enumerate()
                                .map(|(k, a3_coef_ijk)| {
                                    self.x[k]
                                        * a3_coef_ijk
                                        * J3cTerm::new::<2>(
                                            (self.m[i] * self.m[j] * self.m[k]).cbrt(),
                                        )
                                        .t0d0(eta)
                                })
                                .sum::<f64>()
                    })
                    .sum::<f64>()
                    * self.dens2temp3
            })
            .collect();
        zip(a2, a3)
            .map(move |(a2, a3)| coef_a2 * a2 + coef_a3 * a3)
            .collect()
    }
}
/// PdqTerm
const ALPHA: f64 = 1.19374;
#[derive(Clone)]
struct PdqTerm {
    temp: f64,
    dens1temp2: f64, // dens^1 / temp^2
    dens2temp3: f64, // dens^2 / temp^3
    // a2term
    vec_epsilon_ij: Vec<f64>,
    vec_a2_coef: Vec<f64>,
    vec_j2a: Vec<J2aTerm>,
    vec_j2b: Vec<J2bTerm>,
    // a3term
    vec_a3_coef: Vec<f64>,
    vec_j3c: Vec<J3cTerm>,
    // calc mu_k
    x_dd: Vec<f64>,
    x_qq: Vec<f64>,
    m_dd: Vec<f64>,
    m_qq: Vec<f64>,
    epsilon_dd: Vec<f64>, // for mat_epsilon_ij
    epsilon_qq: Vec<f64>, // for mat_epsilon_ij
    mat_a2_coef: Vec<Vec<f64>>,
    mat_a3_coef_ddq: Vec<Vec<Vec<f64>>>,
    mat_a3_coef_dqq: Vec<Vec<Vec<f64>>>,
    a2t0d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t0d4: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t1d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a2t2d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t0d4: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d0: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d1: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d2: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t1d3: (f64, f64), // cached variables -> macro_rules! fn_polar
    a3t2d0: (f64, f64), // cached variables -> macro_rules! fn_polar
}
fn_polar!(PdqTerm);
impl PdqTerm {
    fn new(
        (x_dd, x_qq): (&[f64], &[f64]),
        (m_dd, m_qq): (&[f64], &[f64]),
        (sigma_dd, sigma_qq): (&[f64], &[f64]),
        (epsilon_dd, epsilon_qq): (&[f64], &[f64]),
        (mu, q): (&[f64], &[f64]),
    ) -> Self {
        let mu2_plus: Vec<f64> = mu
            .iter()
            .zip(m_dd.iter())
            .map(|(mu, m)| mu.powi(2) * 1E-19 / m / K)
            .collect();
        let m_dd: Vec<f64> = m_dd.iter().map(|m| m.min(2.0)).collect();
        let q2_plus: Vec<f64> = q
            .iter()
            .zip(m_qq.iter())
            .map(|(q, m)| q.powi(2) * 1E-19 / m / K)
            .collect();
        let m_qq: Vec<f64> = m_qq.iter().map(|m| m.min(2.0)).collect();
        let mat_a2_coef: Vec<Vec<f64>> = mu2_plus
            .iter()
            .zip(sigma_dd)
            .map(|(mu2_i, sigma_i)| {
                q2_plus
                    .iter()
                    .zip(sigma_qq)
                    .map(|(q2_j, sigma_j)| -72.0 * PI * mu2_i * q2_j / (sigma_i + sigma_j).powi(5))
                    .collect()
            })
            .collect();
        let mat_a3_coef_ddq: Vec<Vec<Vec<f64>>> = mu2_plus
            .iter()
            .zip(sigma_dd)
            .map(|(mu2_i, sigma_i)| {
                mu2_plus
                    .iter()
                    .zip(sigma_dd)
                    .map(|(mu2_j, sigma_j)| {
                        q2_plus
                            .iter()
                            .zip(sigma_qq)
                            .map(|(q2_k, sigma_k)| {
                                // this sign is different from publication
                                // but we think the right sign is positive
                                64.0 * (mu2_i * sigma_i) * (mu2_j * sigma_j) * (q2_k / sigma_k)
                                    / ((sigma_i + sigma_j)
                                        * (sigma_i + sigma_k)
                                        * (sigma_j + sigma_k))
                                        .powi(2)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let mat_a3_coef_dqq: Vec<Vec<Vec<f64>>> = mu2_plus
            .iter()
            .zip(sigma_dd)
            .map(|(mu2_i, sigma_i)| {
                q2_plus
                    .iter()
                    .zip(sigma_qq)
                    .map(|(q2_j, sigma_j)| {
                        q2_plus
                            .iter()
                            .zip(sigma_qq)
                            .map(|(q2_k, sigma_k)| {
                                // this sign is different from publication
                                // but we think the right sign is positive
                                64.0 * ALPHA
                                    * (mu2_i * sigma_i)
                                    * (q2_j / sigma_j)
                                    * (q2_k / sigma_k)
                                    / ((sigma_i + sigma_j)
                                        * (sigma_i + sigma_k)
                                        * (sigma_j + sigma_k))
                                        .powi(2)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let mut vec_epsilon_ij: Vec<f64> = Vec::new();
        let mut vec_a2_coef: Vec<f64> = Vec::new();
        let mut vec_j2a: Vec<J2aTerm> = Vec::new();
        let mut vec_j2b: Vec<J2bTerm> = Vec::new();
        let _ = (0..x_dd.len())
            .map(|i| {
                let _ = (0..x_qq.len())
                    .map(|j| {
                        vec_epsilon_ij.push((epsilon_dd[i] * epsilon_qq[j]).sqrt());
                        vec_a2_coef.push(x_dd[i] * x_qq[j] * mat_a2_coef[i][j]);
                        let m_ij = (m_dd[i] * m_qq[j]).sqrt();
                        vec_j2a.push(J2aTerm::new::<3>(m_ij));
                        vec_j2b.push(J2bTerm::new::<3>(m_ij));
                    })
                    .count();
            })
            .count();
        let mut vec_a3_coef: Vec<f64> = Vec::new();
        let mut vec_j3c: Vec<J3cTerm> = Vec::new();
        let _ = (0..x_qq.len())
            .map(|k| {
                let _ = (0..x_dd.len())
                    .map(|i| {
                        let _ = (0..i)
                            .map(|j| {
                                vec_a3_coef.push(
                                    2.0 * x_dd[i] * x_dd[j] * x_qq[k] * mat_a3_coef_ddq[i][j][k],
                                );
                                vec_j3c
                                    .push(J3cTerm::new::<3>((m_dd[i] * m_dd[j] * m_qq[k]).cbrt()));
                            })
                            .count();
                        vec_a3_coef.push(x_dd[i] * x_dd[i] * x_qq[k] * mat_a3_coef_ddq[i][i][k]);
                        vec_j3c.push(J3cTerm::new::<3>((m_dd[i] * m_dd[i] * m_qq[k]).cbrt()));
                    })
                    .count();
            })
            .count();
        let _ = (0..x_dd.len())
            .map(|i| {
                let _ = (0..x_qq.len())
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                vec_a3_coef.push(
                                    2.0 * x_dd[i] * x_qq[j] * x_qq[k] * mat_a3_coef_dqq[i][j][k],
                                );
                                vec_j3c
                                    .push(J3cTerm::new::<3>((m_dd[i] * m_qq[j] * m_qq[k]).cbrt()));
                            })
                            .count();
                        vec_a3_coef.push(x_dd[i] * x_qq[j] * x_qq[j] * mat_a3_coef_dqq[i][j][j]);
                        vec_j3c.push(J3cTerm::new::<3>((m_dd[i] * m_qq[j] * m_qq[j]).cbrt()));
                    })
                    .count();
            })
            .count();
        Self {
            temp: 0.0,
            dens1temp2: 0.0,
            dens2temp3: 0.0,
            // a2term
            vec_epsilon_ij,
            vec_a2_coef,
            vec_j2a,
            vec_j2b,
            // a3term
            vec_a3_coef,
            vec_j3c,
            // calc mu_k
            x_dd: x_dd.to_owned(),
            x_qq: x_qq.to_owned(),
            m_dd: m_dd.to_owned(),
            m_qq: m_qq.to_owned(),
            epsilon_dd: epsilon_dd.to_owned(),
            epsilon_qq: epsilon_qq.to_owned(),
            mat_a2_coef,
            mat_a3_coef_ddq,
            mat_a3_coef_dqq,
            a2t0d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t0d4: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t1d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a2t2d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t0d4: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d1: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d2: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t1d3: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
            a3t2d0: (0.0, 0.0), // cached variables -> macro_rules! fn_polar
        }
    }
    fn new_fracs(&self, (x_dd, x_qq): (&[f64], &[f64])) -> Self {
        let mut vec_a2_coef: Vec<f64> = Vec::new();
        let _ = (0..x_dd.len())
            .map(|i| {
                let _ = (0..x_qq.len())
                    .map(|j| {
                        vec_a2_coef.push(x_dd[i] * x_qq[j] * self.mat_a2_coef[i][j]);
                    })
                    .count();
            })
            .count();
        let mut vec_a3_coef: Vec<f64> = Vec::new();
        let _ = (0..x_qq.len())
            .map(|k| {
                let _ = (0..x_dd.len())
                    .map(|i| {
                        let _ = (0..i)
                            .map(|j| {
                                vec_a3_coef.push(
                                    2.0 * x_dd[i]
                                        * x_dd[j]
                                        * x_qq[k]
                                        * self.mat_a3_coef_ddq[i][j][k],
                                );
                            })
                            .count();
                        vec_a3_coef
                            .push(x_dd[i] * x_dd[i] * x_qq[k] * self.mat_a3_coef_ddq[i][i][k]);
                    })
                    .count();
            })
            .count();
        let _ = (0..x_dd.len())
            .map(|i| {
                let _ = (0..x_qq.len())
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                vec_a3_coef.push(
                                    2.0 * x_dd[i]
                                        * x_qq[j]
                                        * x_qq[k]
                                        * self.mat_a3_coef_dqq[i][j][k],
                                );
                            })
                            .count();
                        vec_a3_coef
                            .push(x_dd[i] * x_qq[j] * x_qq[j] * self.mat_a3_coef_dqq[i][j][j]);
                    })
                    .count();
            })
            .count();
        Self {
            x_dd: x_dd.to_owned(),
            x_qq: x_qq.to_owned(),
            vec_a2_coef,
            vec_a3_coef,
            ..self.clone()
        }
    }
    fn mu_dd(&mut self, temp: f64, rho_num: f64, eta: f64) -> Vec<f64> {
        self.td_flash(temp, rho_num);
        let coef_a2 = self.a2t0d0(eta) * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
            / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a2: Vec<f64> = self
            .mat_a2_coef
            .iter()
            .enumerate()
            .map(|(i, a2_coef_i)| {
                a2_coef_i
                    .iter()
                    .enumerate()
                    .map(|(j, a2_coef_ij)| {
                        self.x_qq[j]
                            * a2_coef_ij
                            * (J2aTerm::new::<3>((self.m_dd[i] * self.m_qq[j]).sqrt()).t0d0(eta)
                                + (self.epsilon_dd[i] * self.epsilon_qq[j]).sqrt() / self.temp
                                    * J2bTerm::new::<3>((self.m_dd[i] * self.m_qq[j]).sqrt())
                                        .t0d0(eta))
                    })
                    .sum::<f64>()
                    * self.dens1temp2
            })
            .collect();
        let coef_a3 = self.a2t0d0(eta).powi(2) / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a3_ddq: Vec<f64> = self
            .mat_a3_coef_ddq
            .iter()
            .enumerate()
            .map(|(i, mat_a3_coef_i)| {
                2.0 * mat_a3_coef_i
                    .iter()
                    .enumerate()
                    .map(|(j, a3_coef_ij)| {
                        self.x_dd[j]
                            * a3_coef_ij
                                .iter()
                                .enumerate()
                                .map(|(k, a3_coef_ijk)| {
                                    self.x_qq[k]
                                        * a3_coef_ijk
                                        * J3cTerm::new::<3>(
                                            (self.m_dd[i] * self.m_dd[j] * self.m_qq[k]).cbrt(),
                                        )
                                        .t0d0(eta)
                                })
                                .sum::<f64>()
                    })
                    .sum::<f64>()
                    * self.dens2temp3
            })
            .collect();
        let a3_dqq: Vec<f64> = self
            .mat_a3_coef_dqq
            .iter()
            .enumerate()
            .map(|(i, mat_a3_coef_i)| {
                mat_a3_coef_i
                    .iter()
                    .enumerate()
                    .map(|(j, a3_coef_ij)| {
                        self.x_qq[j]
                            * a3_coef_ij
                                .iter()
                                .enumerate()
                                .map(|(k, a3_coef_ijk)| {
                                    self.x_qq[k]
                                        * a3_coef_ijk
                                        * J3cTerm::new::<3>(
                                            (self.m_dd[i] * self.m_qq[j] * self.m_qq[k]).cbrt(),
                                        )
                                        .t0d0(eta)
                                })
                                .sum::<f64>()
                    })
                    .sum::<f64>()
                    * self.dens2temp3
            })
            .collect();
        a2.into_iter()
            .zip(a3_ddq)
            .zip(a3_dqq)
            .map(move |((a2, a3_ddq), a3_dqq)| coef_a2 * a2 + coef_a3 * (a3_ddq + a3_dqq))
            .collect()
    }
    fn mu_qq(&mut self, temp: f64, rho_num: f64, eta: f64) -> Vec<f64> {
        self.td_flash(temp, rho_num);
        let coef_a2 = self.a2t0d0(eta) * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
            / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a2: Vec<f64> = (0..self.mat_a2_coef[0].len())
            .map(|j| {
                self.mat_a2_coef
                    .iter()
                    .enumerate()
                    .map(|(i, mat_a2_coef_i)| {
                        self.x_dd[i]
                            * mat_a2_coef_i[j]
                            * (J2aTerm::new::<3>((self.m_dd[i] * self.m_qq[j]).sqrt()).t0d0(eta)
                                + (self.epsilon_dd[i] * self.epsilon_qq[j]).sqrt() / self.temp
                                    * J2bTerm::new::<3>((self.m_dd[i] * self.m_qq[j]).sqrt())
                                        .t0d0(eta))
                    })
                    .sum::<f64>()
                    * self.dens1temp2
            })
            .collect();
        let coef_a3 = self.a2t0d0(eta).powi(2) / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2);
        let a3_ddq: Vec<f64> = (0..self.mat_a3_coef_ddq[0].len())
            .map(|k| {
                self.mat_a3_coef_ddq
                    .iter()
                    .enumerate()
                    .map(|(i, mat_a3_coef_i)| {
                        self.x_dd[i]
                            * mat_a3_coef_i
                                .iter()
                                .enumerate()
                                .map(|(j, mat_a3_coef_ij)| {
                                    self.x_dd[j]
                                        * mat_a3_coef_ij[k]
                                        * J3cTerm::new::<3>(
                                            (self.m_dd[i] * self.m_dd[j] * self.m_qq[k]).cbrt(),
                                        )
                                        .t0d0(eta)
                                })
                                .sum::<f64>()
                    })
                    .sum::<f64>()
                    * self.dens2temp3
            })
            .collect();
        let a3_dqq: Vec<f64> = (0..self.mat_a3_coef_dqq[0].len())
            .map(|k| {
                2.0 * self
                    .mat_a3_coef_dqq
                    .iter()
                    .enumerate()
                    .map(|(i, mat_a3_coef_i)| {
                        self.x_dd[i]
                            * mat_a3_coef_i
                                .iter()
                                .enumerate()
                                .map(|(j, mat_a3_coef_ij)| {
                                    self.x_qq[j]
                                        * mat_a3_coef_ij[k]
                                        * J3cTerm::new::<3>(
                                            (self.m_dd[i] * self.m_qq[j] * self.m_qq[k]).cbrt(),
                                        )
                                        .t0d0(eta)
                                })
                                .sum::<f64>()
                    })
                    .sum::<f64>()
                    * self.dens2temp3
            })
            .collect();
        a2.into_iter()
            .zip(a3_ddq)
            .zip(a3_dqq)
            .map(move |((a2, a3_ddq), a3_dqq)| coef_a2 * a2 + coef_a3 * (a3_ddq + a3_dqq))
            .collect()
    }
}
