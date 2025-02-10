use thiserror::Error;
#[derive(Debug, Error)]
enum RkErr {
    #[error("t_flash diverge")]
    NotConvForT,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
    #[error("property only in double phase")]
    OnlyInDoublePhase,
}
const ZC: f64 = 1.0 / 3.0;
fn ac_coef() -> &'static f64 {
    static AC_COEF: OnceLock<f64> = OnceLock::new();
    AC_COEF.get_or_init(|| 4_f64.cbrt() / (2.0 - 4_f64.cbrt()) / 9.0)
}
fn bc_coef() -> &'static f64 {
    static BC_COEF: OnceLock<f64> = OnceLock::new();
    BC_COEF.get_or_init(|| (2_f64.cbrt() - 1.0) / 3.0)
}
use crate::algorithms::shengjin_roots;
use crate::f64consts::R;
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::sync::OnceLock;
/// Rk EOS
/// ```
/// use thermolib::Rk;
/// let crit_t = 430.64; // critical temperature of sulfur dioxide // K
/// let crit_p = 7886600.0; // critical pressure of sulfur dioxide // Pa
/// let mut fluid = Rk::new_fluid(crit_t, crit_p);
/// fluid.t_flash(273.15).unwrap();
/// assert_eq!(fluid.p_s().unwrap().round(), 282859.0);
/// assert_eq!(fluid.rho_v().unwrap().round(), 130.0);
/// assert_eq!(fluid.rho_l().unwrap().round(), 19422.0);
/// fluid.tp_flash(273.15, 0.1e6);
/// assert_eq!(fluid.rho().unwrap().round(), 45.0);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)]
pub struct Rk {
    T: f64,
    p: f64,
    a: f64,
    b: f64,
    A: f64,
    B: f64,
    Z: f64,
    Zv: f64,
    Zl: f64,
    pc: f64,
    is_single_phase: bool,
}
impl Rk {
    pub fn new_fluid(temp_c: f64, pc: f64) -> Self {
        let mut rk = Rk {
            a: ac_coef() / pc * R.powi(2) * temp_c.powf(2.5),
            b: bc_coef() / pc * R * temp_c,
            T: 0.0,
            p: 0.0,
            A: 0.0,
            B: 0.0,
            Z: 0.0,
            Zv: 0.0,
            Zl: 0.0,
            pc,
            is_single_phase: false,
        };
        rk.tp_flash(273.15, 0.1E6);
        rk
    }
}
impl Rk {
    fn calc_root(&mut self) {
        self.A = self.a * self.p / R.powi(2) / self.T.powf(2.5);
        self.B = self.b * self.p / R / self.T;
        let (zv, zl) = shengjin_roots(-1.0, self.A - self.B - self.B.powi(2), -self.A * self.B);
        if zl.is_sign_negative() {
            self.Z = zv;
            self.is_single_phase = true;
        } else {
            (self.Zv, self.Zl) = (zv, zl);
            self.is_single_phase = false;
        }
    }
    fn calc_lnfp(&self, z: f64) -> f64 {
        z - 1.0 - (z - self.B).ln() + self.A / self.B * (z / (z + self.B)).abs().ln()
    }
    fn calc_diff_lnfpvl(&mut self) -> f64 {
        self.calc_root();
        if self.is_single_phase {
            if self.Z > ZC {
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
impl Rk {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(temp_c: f64, pc: f64) -> Self {
        Self::new_fluid(temp_c, pc)
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
            Err(anyhow!(RkErr::NotConvForT))
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
            Err(anyhow!(RkErr::OnlyInSinglePhase))
        }
    }
    pub fn p(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.p)
        } else {
            Err(anyhow!(RkErr::OnlyInSinglePhase))
        }
    }
    pub fn rho(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.p / (self.Z * R * self.T))
        } else {
            Err(anyhow!(RkErr::OnlyInSinglePhase))
        }
    }
    #[allow(non_snake_case)]
    pub fn T_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::OnlyInDoublePhase))
        } else {
            Ok(self.T)
        }
    }
    pub fn p_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::OnlyInDoublePhase))
        } else {
            Ok(self.p)
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::OnlyInDoublePhase))
        } else {
            Ok(self.p / (self.Zv * R * self.T))
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::OnlyInDoublePhase))
        } else {
            Ok(self.p / (self.Zl * R * self.T))
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_rk() {
        let crit_t = 430.64; // critical temperature of sulfur dioxide // K
        let crit_p = 7886600.0; // critical pressure of sulfur dioxide // Pa
        let mut fluid = Rk::new_fluid(crit_t, crit_p);
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
