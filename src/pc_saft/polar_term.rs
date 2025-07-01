use crate::f64consts::K;
use std::f64::consts::PI;
const PI2: f64 = PI * PI;
pub struct PolarTerm {
    qq: Option<PqqTerm>,
    dd: Option<PddTerm>,
    dq: Option<PdqTerm>,
}
impl PolarTerm {
    pub fn new(x: &[f64], m: &[f64], sigma: &[f64], epsilon: &[f64], p: &[f64], n: &[i32]) -> Self {
        // qq
        let n_q: Vec<usize> = n
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n == 4)
            .map(|(i, _)| i)
            .collect();
        let x_q: Vec<f64> = n_q.iter().map(|&i| x[i]).collect();
        let m_q: Vec<f64> = n_q.iter().map(|&i| m[i]).collect();
        let sigma_q: Vec<f64> = n_q.iter().map(|&i| sigma[i]).collect();
        let epsilon_q: Vec<f64> = n_q.iter().map(|&i| epsilon[i]).collect();
        let p_q: Vec<f64> = n_q.iter().map(|&i| p[i]).collect();
        let qq = if n_q.is_empty() {
            None
        } else {
            Some(PqqTerm::new(&x_q, &m_q, &sigma_q, &epsilon_q, &p_q))
        };
        // dd
        let n_d: Vec<usize> = n
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n == 2)
            .map(|(i, _)| i)
            .collect();
        let x_d: Vec<f64> = n_d.iter().map(|&i| x[i]).collect();
        let m_d: Vec<f64> = n_d.iter().map(|&i| m[i]).collect();
        let sigma_d: Vec<f64> = n_d.iter().map(|&i| sigma[i]).collect();
        let epsilon_d: Vec<f64> = n_d.iter().map(|&i| epsilon[i]).collect();
        let p_d: Vec<f64> = n_d.iter().map(|&i| p[i]).collect();
        let dd = if n_d.is_empty() {
            None
        } else {
            Some(PddTerm::new(&x_d, &m_d, &sigma_d, &epsilon_d, &p_d))
        };
        // dq
        Self {
            qq,
            dd,
            dq: if n_d.is_empty() || n_q.is_empty() {
                None
            } else {
                Some(PdqTerm::new(
                    &x_q, &m_q, &sigma_q, &epsilon_q, &p_q, &x_d, &m_d, &sigma_d, &epsilon_d, &p_d,
                ))
            },
        }
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
            fn td_flash(&mut self, temp: f64, rho_num: f64) {
                self.temp = temp;
                self.dens1temp2 = rho_num / temp.powi(2);
                self.dens2temp3 = rho_num.powi(2) / temp.powi(3);
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
            fn a2t0d0(&mut self, eta: f64) -> f64 {
                if eta != self.a2t0d0.0 {
                    self.a2t0d0 = (
                        eta,
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.j2a.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
                        self.a2_coef
                            .iter()
                            .zip(self.epsilon_ij.iter())
                            .zip(self.j2a.iter_mut())
                            .zip(self.j2b.iter_mut())
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
            fn a3t0d0(&mut self, eta: f64) -> f64 {
                if eta != self.a3t0d0.0 {
                    self.a3t0d0 = (
                        eta,
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
                        self.a3_coef
                            .iter()
                            .zip(self.j3c.iter_mut())
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
const DQ_A00: f64 = 0.6970590;
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
struct PqqTerm {
    temp: f64,
    dens1temp2: f64, // dens^1 / temp^2
    dens2temp3: f64, // dens^2 / temp^3
    epsilon_ij: Vec<f64>,
    a2_coef: Vec<f64>,
    a3_coef: Vec<f64>,
    j2a: Vec<J2aTerm>,
    j2b: Vec<J2bTerm>,
    j3c: Vec<J3cTerm>,
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
        let n = x.len();
        let q2_plus: Vec<f64> = x
            .iter()
            .zip(m.iter())
            .zip(q.iter())
            .map(|((x, m), q)| x / m / K * q.powi(2) * 1E-19)
            .collect();
        let mut epsilon_ij: Vec<f64> = Vec::new();
        let mut a2_coef: Vec<f64> = Vec::new();
        let mut j2a: Vec<J2aTerm> = Vec::new();
        let mut j2b: Vec<J2bTerm> = Vec::new();
        let _ = (0..n)
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        epsilon_ij.push((epsilon[i] * epsilon[j]).sqrt());
                        a2_coef.push(
                            -72.0 * PI * q2_plus[i] * q2_plus[j] / (sigma[i] + sigma[j]).powi(7),
                        );
                        let m_ij = (m[i] * m[j]).sqrt().min(2.0);
                        j2a.push(J2aTerm::new::<4>(m_ij));
                        j2b.push(J2bTerm::new::<4>(m_ij));
                    })
                    .count();
                epsilon_ij.push(epsilon[i]);
                a2_coef.push(-0.5625 * PI * q2_plus[i].powi(2) / sigma[i].powi(7));
                j2a.push(J2aTerm::new::<4>(m[i].min(2.0)));
                j2b.push(J2bTerm::new::<4>(m[i].min(2.0)));
            })
            .count();
        let mut a3_coef: Vec<f64> = Vec::new();
        let mut j3c: Vec<J3cTerm> = Vec::new();
        let _ = (0..n)
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                a3_coef.push(
                                    1728.0 * PI2 * q2_plus[i] * q2_plus[j] * q2_plus[k]
                                        / ((sigma[i] + sigma[j])
                                            * (sigma[i] + sigma[k])
                                            * (sigma[j] + sigma[k]))
                                            .powi(3),
                                );
                                j3c.push(J3cTerm::new::<4>((m[i] * m[j] * m[k]).cbrt().min(2.0)));
                            })
                            .count();
                        a3_coef.push(
                            216.0 * PI2 * (q2_plus[i] * q2_plus[j].powi(2))
                                / ((sigma[i] + sigma[j]).powi(6) * sigma[j].powi(3)),
                        );
                        j3c.push(J3cTerm::new::<4>((m[i] * m[j] * m[j]).cbrt().min(2.0)));
                    })
                    .count();
                a3_coef.push(0.5625 * PI2 * q2_plus[i].powi(3) / sigma[i].powi(9));
                j3c.push(J3cTerm::new::<4>(m[i].min(2.0)));
            })
            .count();
        Self {
            temp: 0.0,
            dens1temp2: 0.0,
            dens2temp3: 0.0,
            epsilon_ij,
            a2_coef,
            a3_coef,
            j2a,
            j2b,
            j3c,
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
}
/// PddTerm
struct PddTerm {
    temp: f64,
    dens1temp2: f64, // dens^1 / temp^2
    dens2temp3: f64, // dens^2 / temp^3
    epsilon_ij: Vec<f64>,
    a2_coef: Vec<f64>,
    a3_coef: Vec<f64>,
    j2a: Vec<J2aTerm>,
    j2b: Vec<J2bTerm>,
    j3c: Vec<J3cTerm>,
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
    pub fn new(x: &[f64], m: &[f64], sigma: &[f64], epsilon: &[f64], u: &[f64]) -> Self {
        let n = x.len();
        let u2_plus: Vec<f64> = x
            .iter()
            .zip(m.iter())
            .zip(u.iter())
            .map(|((x, m), u)| x / m / K * u.powi(2) * 1E-19)
            .collect();
        let mut epsilon_ij: Vec<f64> = Vec::new();
        let mut a2_coef: Vec<f64> = Vec::new();
        let mut j2a: Vec<J2aTerm> = Vec::new();
        let mut j2b: Vec<J2bTerm> = Vec::new();
        let _ = (0..n)
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        epsilon_ij.push((epsilon[i] * epsilon[j]).sqrt());
                        a2_coef.push(
                            -16.0 * PI * u2_plus[i] * u2_plus[j] / (sigma[i] + sigma[j]).powi(3),
                        );
                        let m_ij = (m[i] * m[j]).sqrt().min(2.0);
                        j2a.push(J2aTerm::new::<2>(m_ij));
                        j2b.push(J2bTerm::new::<2>(m_ij));
                    })
                    .count();
                epsilon_ij.push(epsilon[i]);
                a2_coef.push(-PI * u2_plus[i].powi(2) / sigma[i].powi(3));
                j2a.push(J2aTerm::new::<2>(m[i].min(2.0)));
                j2b.push(J2bTerm::new::<2>(m[i].min(2.0)));
            })
            .count();
        let mut a3_coef: Vec<f64> = Vec::new();
        let mut j3c: Vec<J3cTerm> = Vec::new();
        let _ = (0..n)
            .map(|i| {
                let _ = (0..i)
                    .map(|j| {
                        let _ = (0..j)
                            .map(|k| {
                                a3_coef.push(
                                    -64.0 * PI2 * u2_plus[i] * u2_plus[j] * u2_plus[k]
                                        / (sigma[i] + sigma[j])
                                        / (sigma[i] + sigma[k])
                                        / (sigma[j] + sigma[k]),
                                );
                                j3c.push(J3cTerm::new::<2>((m[i] * m[j] * m[k]).cbrt().min(2.0)));
                            })
                            .count();
                        a3_coef.push(
                            -32.0 * PI2 * u2_plus[i] * u2_plus[j].powi(2)
                                / ((sigma[i] + sigma[j]).powi(2) * sigma[j]),
                        );
                        j3c.push(J3cTerm::new::<2>((m[i] * m[j] * m[j]).cbrt().min(2.0)));
                    })
                    .count();
                a3_coef.push(-4.0 / 3.0 * PI2 * u2_plus[i].powi(3) / sigma[i].powi(3));
                j3c.push(J3cTerm::new::<2>(m[i].min(2.0)));
            })
            .count();
        Self {
            temp: 0.0,
            dens1temp2: 0.0,
            dens2temp3: 0.0,
            epsilon_ij,
            a2_coef,
            a3_coef,
            j2a,
            j2b,
            j3c,
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
}
/// PdqTerm
const ALPHA: f64 = 1.19374;
struct PdqTerm {
    temp: f64,
    dens1temp2: f64, // dens^1 / temp^2
    dens2temp3: f64, // dens^2 / temp^3
    epsilon_ij: Vec<f64>,
    a2_coef: Vec<f64>,
    a3_coef: Vec<f64>,
    j2a: Vec<J2aTerm>,
    j2b: Vec<J2bTerm>,
    j3c: Vec<J3cTerm>,
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
    #[allow(clippy::too_many_arguments)]
    fn new(
        x_q: &[f64],
        m_q: &[f64],
        sigma_q: &[f64],
        epsilon_q: &[f64],
        p_q: &[f64],
        x_d: &[f64],
        m_d: &[f64],
        sigma_d: &[f64],
        epsilon_d: &[f64],
        p_d: &[f64],
    ) -> Self {
        let q2_plus: Vec<f64> = x_q
            .iter()
            .zip(m_q.iter())
            .zip(p_q.iter())
            .map(|((x, m), q)| x / m / K * q.powi(2) * 1E-19)
            .collect();
        let q2_len = q2_plus.len();
        let u2_plus: Vec<f64> = x_d
            .iter()
            .zip(m_d.iter())
            .zip(p_d.iter())
            .map(|((x, m), u)| x / m / K * u.powi(2) * 1E-19)
            .collect();
        let u2_len = u2_plus.len();
        let mut epsilon_ij: Vec<f64> = Vec::new();
        let mut a2_coef: Vec<f64> = Vec::new();
        let mut j2a: Vec<J2aTerm> = Vec::new();
        let mut j2b: Vec<J2bTerm> = Vec::new();
        let _ = (0..u2_len).map(|i| {
            let _ = (0..q2_len).map(|j| {
                epsilon_ij.push((epsilon_d[i] * epsilon_q[j]).sqrt());
                a2_coef
                    .push(-72.0 * PI * u2_plus[i] * q2_plus[j] / (sigma_d[i] + sigma_q[j]).powi(5));
                let m_ij = (m_d[i] * m_q[j]).sqrt().min(2.0);
                j2a.push(J2aTerm::new::<3>(m_ij));
                j2b.push(J2bTerm::new::<3>(m_ij));
            });
        });
        let mut a3_coef: Vec<f64> = Vec::new();
        let mut j3c: Vec<J3cTerm> = Vec::new();
        let _ = (0..q2_len)
            .map(|k| {
                let _ = (0..u2_len)
                    .map(|i| {
                        let _ = (0..i)
                            .map(|j| {
                                a3_coef.push(
                                    -128.0
                                        * (u2_plus[i] * sigma_d[i])
                                        * (u2_plus[j] * sigma_d[j])
                                        * (q2_plus[k] / sigma_q[k])
                                        / ((sigma_d[i] + sigma_d[j])
                                            * (sigma_d[i] + sigma_q[k])
                                            * (sigma_d[j] + sigma_q[k]))
                                            .powi(2),
                                );
                                j3c.push(J3cTerm::new::<3>(
                                    (m_d[i] * m_d[j] * m_q[k]).cbrt().min(2.0),
                                ));
                            })
                            .count();
                        a3_coef.push(
                            -16.0 * (u2_plus[i] * sigma_d[i]).powi(2) * (q2_plus[k] / sigma_q[k])
                                / sigma_d[i].powi(2)
                                / (sigma_d[i] + sigma_q[k]).powi(4),
                        );
                        j3c.push(J3cTerm::new::<3>(
                            (m_d[i] * m_d[i] * m_q[k]).cbrt().min(2.0),
                        ));
                    })
                    .count();
            })
            .count();
        let _ = (0..u2_len).map(|i| {
            let _ = (0..q2_len).map(|j| {
                let _ = (0..j).map(|k| {
                    a3_coef.push(
                        -128.0
                            * ALPHA
                            * (u2_plus[i] * sigma_d[i])
                            * (q2_plus[j] / sigma_q[j])
                            * (q2_plus[k] / sigma_q[k])
                            / ((sigma_d[i] + sigma_q[j])
                                * (sigma_d[i] + sigma_q[k])
                                * (sigma_q[j] + sigma_q[k]))
                                .powi(2),
                    );
                    j3c.push(J3cTerm::new::<3>(
                        (m_d[i] * m_q[j] * m_q[k]).cbrt().min(2.0),
                    ));
                });
                a3_coef.push(
                    -16.0 * ALPHA * (u2_plus[i] * sigma_d[i]) * (q2_plus[j] / sigma_q[j]).powi(2)
                        / (sigma_d[i] + sigma_q[i]).powi(4)
                        / sigma_q[j].powi(2),
                );
                j3c.push(J3cTerm::new::<3>(
                    (m_d[i] * m_q[j] * m_q[j]).cbrt().min(2.0),
                ));
            });
        });
        Self {
            temp: 0.0,
            dens1temp2: 0.0,
            dens2temp3: 0.0,
            epsilon_ij,
            a2_coef,
            a3_coef,
            j2a,
            j2b,
            j3c,
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
}
