use thiserror::Error;
#[derive(Debug, Error)]
enum PrErr {
    #[error("t_flash diverge")]
    NotConvForT,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
    #[error("property only in double phase")]
    OnlyInDoublePhase,
}
use std::f64::consts::SQRT_2;
const SQRT2_ADD_POS1: f64 = SQRT_2 + 1.0;
const SQRT2_ADD_NEG1: f64 = SQRT_2 - 1.0;
const R: f64 = 8.314462618;
fn c() -> &'static f64 {
    static C: OnceLock<f64> = OnceLock::new();
    C.get_or_init(|| ((6.0 * SQRT_2 + 8.0).cbrt() - (6.0 * SQRT_2 - 8.0).cbrt() - 1.0) / 3.0)
}
fn zc() -> &'static f64 {
    static AC_COEF: OnceLock<f64> = OnceLock::new();
    AC_COEF
        .get_or_init(|| (1.0 - 2.0 * c() - c().powi(2)) / (1.0 + c()) / (1.0 - c()).powi(2) / 2.0)
}
fn ac_coef() -> &'static f64 {
    static AC_COEF: OnceLock<f64> = OnceLock::new();
    AC_COEF.get_or_init(|| {
        ((1.0 + 2.0 * c() - c().powi(2)) / (1.0 + c()) / (1.0 - c()).powi(2) / 2.0).powi(2)
            * (1.0 - 2.0 * c() - c().powi(2))
    })
}
fn bc_coef() -> &'static f64 {
    static BC_COEF: OnceLock<f64> = OnceLock::new();
    BC_COEF.get_or_init(|| {
        c() * (1.0 - 2.0 * c() - c().powi(2)) / (1.0 + c()) / (1.0 - c()).powi(2) / 2.0
    })
}
use crate::algorithms::shengjin_roots;
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::sync::OnceLock;
/// Pr EOS
/// ```
/// use thermolib::Pr;
/// let crit_t = 430.64; // critical temperature of sulfur dioxide // K
/// let crit_p = 7886600.0; // critical pressure of sulfur dioxide // Pa
/// let acentric_factor = 0.256;
/// let mut fluid = Pr::new_fluid(crit_t, crit_p, acentric_factor);
/// fluid.t_flash(273.15).unwrap();
/// assert_eq!(fluid.p_s().unwrap().round(), 156174.0);
/// assert_eq!(fluid.rho_v().unwrap().round(), 71.0);
/// assert_eq!(fluid.rho_l().unwrap().round(), 22695.0);
/// fluid.tp_flash(273.15, 0.1e6);
/// assert_eq!(fluid.rho().unwrap().round(), 45.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)]
pub struct Pr {
    kappa: f64,
    ac: f64,
    bc: f64,
    T: f64,
    p: f64,
    A: f64,
    B: f64,
    Z: f64,
    Zv: f64,
    Zl: f64,
    Tc: f64,
    pc: f64,
    is_single_phase: bool,
}
impl Pr {
    pub fn new_fluid(temp_c: f64, pc: f64, omega: f64) -> Self {
        let mut pr = Pr {
            kappa: 0.37464 + 1.54226 * omega - 0.26992 * omega.powi(2),
            ac: ac_coef() / pc * (R * temp_c).powi(2),
            bc: bc_coef() / pc * R * temp_c,
            T: 0.0,
            p: 0.0,
            A: 0.0,
            B: 0.0,
            Z: 0.0,
            Zv: 0.0,
            Zl: 0.0,
            Tc: temp_c,
            pc,
            is_single_phase: false,
        };
        pr.tp_flash(273.15, 0.1E6);
        pr
    }
}
impl Pr {
    fn calc_root(&mut self) {
        self.A = self.ac * (1.0 + self.kappa * (1.0 - (self.T / self.Tc).sqrt())).powi(2) * self.p
            / (R * self.T).powi(2);
        self.B = self.bc * self.p / (R * self.T);
        let (zv, zl) = shengjin_roots(
            self.B - 1.0,
            self.A - 3.0 * self.B.powi(2) - 2.0 * self.B,
            -self.A * self.B + self.B.powi(2) + self.B.powi(3),
        );
        if zl == 0.0 {
            self.Z = zv;
            self.is_single_phase = true;
        } else {
            (self.Zv, self.Zl) = (zv, zl);
            self.is_single_phase = false;
        }
    }
    fn calc_lnfp(&self, z: f64) -> f64 {
        z - 1.0
            - (z - self.B).ln()
            - self.A / (2.0 * SQRT_2 * self.B)
                * ((z + SQRT2_ADD_POS1 * self.B) / (z - SQRT2_ADD_NEG1 * self.B))
                    .abs()
                    .ln()
    }
    fn calc_diff_lnfpvl(&mut self) -> f64 {
        self.calc_root();
        if self.is_single_phase {
            if self.Z > *zc() {
                self.calc_lnfp(self.Z)
            } else {
                -self.calc_lnfp(self.Z)
            }
        } else {
            self.calc_lnfp(self.Zv) - self.calc_lnfp(self.Zl)
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
impl Pr {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(temp_c: f64, pc: f64, omega: f64) -> Self {
        Self::new_fluid(temp_c, pc, omega)
    }
    pub fn t_flash(&mut self, temp: f64) -> anyhow::Result<()> {
        self.T = temp;
        let ps_max = self.pc - 1.0;
        let ps = crate::algorithms::brent_zero(
            |ps| {
                self.p = ps;
                self.calc_diff_lnfpvl()
            },
            1.0,
            ps_max,
        );
        if ps.is_nan() {
            Err(anyhow!(PrErr::NotConvForT))
        } else {
            Ok(())
        }
    }
    pub fn tp_flash(&mut self, temp: f64, p: f64) {
        self.T = temp;
        self.p = p;
        self.calc_root();
        if !self.is_single_phase {
            let lnfpv = self.calc_lnfp(self.Zv);
            let lnfpl = self.calc_lnfp(self.Zl);
            if lnfpv < lnfpl {
                self.Z = self.Zv;
            } else {
                self.Z = self.Zl;
            }
            self.is_single_phase = true;
        }
    }
    #[allow(non_snake_case)]
    pub fn T(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.T)
        } else {
            Err(anyhow!(PrErr::OnlyInSinglePhase))
        }
    }
    pub fn p(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.p)
        } else {
            Err(anyhow!(PrErr::OnlyInSinglePhase))
        }
    }
    pub fn rho(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.p / (self.Z * R * self.T))
        } else {
            Err(anyhow!(PrErr::OnlyInSinglePhase))
        }
    }
    #[allow(non_snake_case)]
    pub fn T_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PrErr::OnlyInDoublePhase))
        } else {
            Ok(self.T)
        }
    }
    pub fn p_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PrErr::OnlyInDoublePhase))
        } else {
            Ok(self.p)
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PrErr::OnlyInDoublePhase))
        } else {
            Ok(self.p / (self.Zv * R * self.T))
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(PrErr::OnlyInDoublePhase))
        } else {
            Ok(self.p / (self.Zl * R * self.T))
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pr() {
        let crit_t = 430.64; // critical temperature of sulfur dioxide // K
        let crit_p = 7886600.0; // critical pressure of sulfur dioxide // Pa
        let acentric_factor = 0.256;
        let mut fluid = Pr::new_fluid(crit_t, crit_p, acentric_factor);
        let temp_min = (0.6 * crit_t).floor() as u32;
        let temp_max = crit_t.ceil() as u32;
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
    }
}
