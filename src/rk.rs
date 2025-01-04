use thiserror::Error;
#[derive(Debug, Error)]
enum RkErr {
    #[error("t_flash diverge")]
    NotConvForT,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
    #[error("property not in single phase")]
    NotInSinglePhase,
}
use crate::algorithms::shengjin_roots;
use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
/// Rk EOS
/// ```
/// use thermolib::Rk;
/// let Tc = 430.64; // K
/// let pc = 7886600.0; // Pa
/// let M = 0.064064; // kg/mol
/// let mut SO2 = Rk::new_fluid(Tc, pc, M);
/// let _ = SO2.set_molar_unit();
/// if let Ok(_) = SO2.t_flash(273.15) {
///     println!("T_s={}", SO2.T_s().unwrap());
///     println!("p_s={}", SO2.p_s().unwrap());
///     println!("rho_v={}", SO2.rho_v().unwrap());
///     println!("rho_l={}", SO2.rho_l().unwrap());
/// }
/// SO2.tp_flash(273.15, 0.1e6);
/// println!("T={}", SO2.T().unwrap());
/// println!("p={}", SO2.p().unwrap());
/// println!("rho={}", SO2.rho().unwrap());
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)]
pub struct Rk {
    Tc: f64,
    pc: f64,
    T: f64,
    p: f64,
    a: f64,
    b: f64,
    A: f64,
    B: f64,
    R: f64,
    M: f64,
    Z: f64,
    Zc: f64,
    Zv: f64,
    Zl: f64,
    is_single_phase: bool,
}
impl Rk {
    #[allow(non_snake_case)]
    pub fn new_fluid(Tc: f64, pc: f64, M: f64) -> Self {
        let R = 8.314462618;
        let Zc = 1.0 / 3.0;
        let mut rk = Rk {
            Zc,
            Tc,
            pc,
            R,
            M,
            a: 0.42748 * R.powi(2) * Tc.powf(2.5) / pc,
            b: 0.08664 * R * Tc / pc,
            T: 0.0,
            p: 0.0,
            A: 0.0,
            B: 0.0,
            Z: 0.0,
            Zv: 0.0,
            Zl: 0.0,
            is_single_phase: false,
        };
        rk.tp_flash(273.15, 0.1E6);
        rk
    }
}
#[allow(non_snake_case)]
impl Rk {
    fn calc_root(&mut self) {
        self.A = self.a * self.p / self.R.powi(2) / self.T.powf(2.5);
        self.B = self.b * self.p / self.R / self.T;
        let (Zv, Zl) = shengjin_roots(-1.0, self.A - self.B - self.B.powi(2), -self.A * self.B);
        if Zl == 0.0 {
            self.Z = Zv;
            self.is_single_phase = true;
        } else {
            (self.Zv, self.Zl) = (Zv, Zl);
            self.is_single_phase = false;
        }
    }
    fn calc_lnfp(&self, Z: f64) -> f64 {
        Z - 1.0 - (Z - self.B).ln() + self.A / self.B * (Z / (Z + self.B)).abs().ln()
    }
    fn calc_diff_lnfpvl(&mut self) -> f64 {
        self.calc_root();
        if self.is_single_phase {
            if self.Z > self.Zc {
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
#[allow(non_snake_case)]
impl Rk {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(Tc: f64, pc: f64, M: f64) -> Self {
        Self::new_fluid(Tc, pc, M)
    }
    pub fn set_molar_unit(&mut self) {
        if self.R > 10.0 {
            self.R *= self.M;
        }
        self.a = 0.42748 * self.R.powi(2) * self.Tc.powf(2.5) / self.pc;
        self.b = 0.08664 * self.R * self.Tc / self.pc;
    }
    pub fn set_mass_unit(&mut self) {
        if self.R < 10.0 {
            self.R /= self.M;
        }
        self.a = 0.42748 * self.R.powi(2) * self.Tc.powf(2.5) / self.pc;
        self.b = 0.08664 * self.R * self.Tc / self.pc;
    }
    pub fn t_flash(&mut self, T: f64) -> anyhow::Result<()> {
        self.T = T;
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
    pub fn tp_flash(&mut self, T: f64, p: f64) {
        self.T = T;
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
            Ok(self.p / (self.Z * self.R * self.T))
        } else {
            Err(anyhow!(RkErr::OnlyInSinglePhase))
        }
    }
    pub fn T_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::NotInSinglePhase))
        } else {
            Ok(self.T)
        }
    }
    pub fn p_s(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::NotInSinglePhase))
        } else {
            Ok(self.p)
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::NotInSinglePhase))
        } else {
            Ok(self.p / (self.Zv * self.R * self.T))
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(RkErr::NotInSinglePhase))
        } else {
            Ok(self.p / (self.Zl * self.R * self.T))
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[allow(non_snake_case)]
    fn test_rk() {
        let Tc: f64 = 430.64; // K
        let pc = 7886600.0; // Pa
        let M = 0.064064; // kg/mol
        let mut SO2 = Rk::new_fluid(Tc, pc, M);
        let Tmin = (0.7 * Tc).floor() as i32;
        let Tmax = Tc.ceil() as i32;
        for T in Tmin..Tmax {
            if let Err(_) = SO2.t_flash(T as f64) {
                panic!();
            }
        }
    }
}
