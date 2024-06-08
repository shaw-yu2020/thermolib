use anyhow::anyhow;
use thiserror::Error;
#[derive(Debug, Error)]
enum SrkErr {
    #[error("T is too min")]
    TisTooMin,
    #[error("T is too max")]
    TisTooMax,
    #[error("t_flash diverge")]
    NotConvForT,
    #[error("property only in single phase")]
    OnlyInSinglePhase,
    #[error("property not in single phase")]
    NotInSinglePhase,
}
use pyo3::{pyclass, pymethods};
/// Srk EOS
#[pyclass]
#[allow(non_snake_case)]
pub struct Srk {
    Tc: f64,
    pc: f64,
    T: f64,
    p: f64,
    m: f64,
    ac: f64,
    bc: f64,
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
#[allow(non_snake_case)]
impl Srk {
    fn calc_root(&mut self) {
        self.A = self.ac * (1.0 + self.m * (1.0 - (self.T / self.Tc).sqrt())).powi(2) * self.p
            / (self.R * self.T).powi(2);
        self.B = self.bc * self.p / (self.R * self.T);
        let c = self.A - self.B - self.B.powi(2);
        let d = -self.A * self.B;
        let A = 1.0 - 3.0 * c;
        let B = -c - 9.0 * d;
        let C = c.powi(2) + 3.0 * d;
        let Delta = B.powi(2) - 4.0 * A * C;
        if Delta.is_sign_negative() {
            let theta3 = ((-2.0 * A - 3.0 * B) / (2.0 * A * A.sqrt())).acos() / 3.0;
            let x1 = (1.0 - 2.0 * A.sqrt() * theta3.cos()) / 3.0;
            let x2 = (1.0 + A.sqrt() * (theta3.cos() + 3_f64.sqrt() * theta3.sin())) / 3.0;
            let x3 = (1.0 + A.sqrt() * (theta3.cos() - 3_f64.sqrt() * theta3.sin())) / 3.0;
            let Zv = x1.max(x2).max(x3);
            let Zl = x1.min(x2).min(x3);
            if Zl.is_sign_negative() {
                self.Z = Zv;
                self.is_single_phase = true;
            } else {
                self.Zv = Zv;
                self.Zl = Zl;
                self.is_single_phase = false;
            }
        } else {
            let Y1 = -A + 1.5 * (-B + Delta.sqrt());
            let Y2 = -A + 1.5 * (-B - Delta.sqrt());
            self.Z = (1.0 - (Y1.cbrt() + Y2.cbrt())) / 3.0;
            self.is_single_phase = true;
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
#[pymethods]
#[allow(non_snake_case)]
impl Srk {
    #[new]
    pub fn new_fluid(Tc: f64, pc: f64, omega: f64, M: f64) -> Self {
        let R = 8.314462618;
        let Zc = 1.0 / 3.0;
        let mut srk = Srk {
            Zc,
            Tc,
            pc,
            R,
            M,
            m: 0.48 + 1.574 * omega - 0.176 * omega.powi(2),
            ac: 0.42747 * (R * Tc).powi(2) / pc,
            bc: 0.08664 * R * Tc / pc,
            T: 0.0,
            p: 0.0,
            A: 0.0,
            B: 0.0,
            Z: 0.0,
            Zv: 0.0,
            Zl: 0.0,
            is_single_phase: false,
        };
        srk.tp_flash(273.15, 0.1E6);
        srk
    }
    pub fn set_mole_unit(&mut self) {
        if self.R > 10.0 {
            self.R *= self.M;
        }
        self.ac = 0.42747 * (self.R * self.Tc).powi(2) / self.pc;
        self.bc = 0.08664 * self.R * self.Tc / self.pc;
    }
    pub fn set_mass_unit(&mut self) {
        if self.R < 10.0 {
            self.R /= self.M;
        }
        self.ac = 0.42747 * (self.R * self.Tc).powi(2) / self.pc;
        self.bc = 0.08664 * self.R * self.Tc / self.pc;
    }
    pub fn t_flash(&mut self, T: f64) -> anyhow::Result<()> {
        self.T = T;
        let mut ps_min = 1.0;
        self.p = ps_min;
        let mut lnfp_pmin = self.calc_diff_lnfpvl();
        let mut ps_max = self.pc - 1.0;
        self.p = ps_max;
        let mut lnfp_pmax = self.calc_diff_lnfpvl();
        if !(lnfp_pmin * lnfp_pmax).is_sign_negative() {
            if T < 0.9 * self.Tc {
                return Err(anyhow!(SrkErr::TisTooMin));
            } else {
                return Err(anyhow!(SrkErr::TisTooMax));
            }
        }
        let mut ps;
        let mut lnfp;
        let mut counter = 0;
        while counter <= 1000 {
            ps = (ps_min + ps_max) / 2.0;
            self.p = ps;
            lnfp = self.calc_diff_lnfpvl();
            if lnfp.abs() < 1E-6 {
                return Ok(());
            } else if (lnfp * lnfp_pmin).is_sign_negative() {
                ps_max = ps;
                lnfp_pmax = lnfp;
            } else if (lnfp * lnfp_pmax).is_sign_negative() {
                ps_min = ps;
                lnfp_pmin = lnfp;
            }
            counter += 1;
        }
        Err(anyhow!(SrkErr::NotConvForT))
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
    pub fn T(&self) -> f64 {
        self.T
    }
    pub fn p(&self) -> f64 {
        self.p
    }
    pub fn rho(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Ok(self.p / (self.Z * self.R * self.T))
        } else {
            Err(anyhow!(SrkErr::OnlyInSinglePhase))
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(SrkErr::NotInSinglePhase))
        } else {
            Ok(self.p / (self.Zv * self.R * self.T))
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        if self.is_single_phase {
            Err(anyhow!(SrkErr::NotInSinglePhase))
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
    fn test_srk() {
        let Tc: f64 = 430.64; // K
        let pc = 7886600.0; // Pa
        let omega = 0.256;
        let M = 0.064064; // kg/mol
        let mut SO2 = Srk::new_fluid(Tc, pc, omega, M);
        let Tmin = (0.7 * Tc).floor() as i32;
        let Tmax = Tc.ceil() as i32;
        for T in Tmin..Tmax {
            if let Err(_) = SO2.t_flash(T as f64) {
                // println!("test_srk panic at {}K", T);
                panic!();
            } else {
                /*
                println!(
                    "test_srk t_flash() at {}K p_s={} rho_v={} rho_l={}",
                    SO2.T(),
                    SO2.p(),
                    SO2.rho_v().unwrap(),
                    SO2.rho_l().unwrap(),
                );
                 */
            }
        }
    }
}
