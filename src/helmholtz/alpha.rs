use super::rho_pr;
use super::{HelmholtzErr, Phase};
use super::{IdealHelmholtz, ResidualHelmholtz};
use super::{PsEqn, RholEqn, RhovEqn};
use anyhow::anyhow;
use pyo3::{pyclass, pymethods};
use serde::Deserialize;
#[pyclass]
#[derive(Deserialize)]
#[allow(non_snake_case)]
pub struct Helmholtz {
    alphai: IdealHelmholtz,
    alphar: ResidualHelmholtz,
    ps: PsEqn,
    rhovs: RhovEqn,
    rhols: RholEqn,
    R: f64,
    M: f64,
    Tc: f64,
    pc: f64,
    rhoc: f64,
    omega: f64,
    #[serde(skip, default)]
    phase: Phase,
}
#[allow(non_snake_case)]
impl Helmholtz {
    fn calc_p(&self, T: f64, rho: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = rho / self.rhoc;
        rho * self.R * T * (1.0 + self.alphar.tau0delta1(tau, delta))
    }
    fn calc_Dp_Drho_T(&self, T: f64, rho: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = rho / self.rhoc;
        (self.R * T)
            * (1.0 + 2.0 * self.alphar.tau0delta1(tau, delta) + self.alphar.tau0delta2(tau, delta))
    }
    fn calc_cv(&self, T: f64, rho: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = rho / self.rhoc;
        self.R * (-self.alphai.tau2(tau, self.Tc) - self.alphar.tau2delta0(tau, delta))
    }
    fn calc_cp(&self, T: f64, rho: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = rho / self.rhoc;
        self.R
            * (-self.alphai.tau2(tau, self.Tc) - self.alphar.tau2delta0(tau, delta)
                + (1.0 + self.alphar.tau0delta1(tau, delta) - self.alphar.tau1delta1(tau, delta))
                    .powi(2)
                    / (1.0
                        + 2.0 * self.alphar.tau0delta1(tau, delta)
                        + self.alphar.tau0delta2(tau, delta)))
    }
    fn calc_w(&self, T: f64, rho: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = rho / self.rhoc;
        let w2 =
            1.0 + 2.0 * self.alphar.tau0delta1(tau, delta) + self.alphar.tau0delta2(tau, delta)
                - (1.0 + self.alphar.tau0delta1(tau, delta) - self.alphar.tau1delta1(tau, delta))
                    .powi(2)
                    / (self.alphai.tau2(tau, self.Tc) + self.alphar.tau2delta0(tau, delta));
        if self.R < 10.0 {
            (self.R / self.M * T * w2).sqrt()
        } else {
            (self.R * T * w2).sqrt()
        }
    }
    fn calc_h(&self, T: f64, rho: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = rho / self.rhoc;
        (self.R * T)
            * (1.0
                + self.alphai.tau1(tau, self.Tc)
                + self.alphar.tau1delta0(tau, delta)
                + self.alphar.tau0delta1(tau, delta))
    }
    fn calc_s(&self, T: f64, rho: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = rho / self.rhoc;
        self.R
            * (self.alphai.tau1(tau, self.Tc) + self.alphar.tau1delta0(tau, delta)
                - self.alphai.tau0(tau, self.Tc)
                - delta.ln()
                - self.alphar.tau0delta0(tau, delta))
    }
    fn J(&self, tau: f64, delta: f64) -> f64 {
        delta * (1.0 + self.alphar.tau0delta1(tau, delta))
    }
    fn K(&self, tau: f64, delta: f64) -> f64 {
        self.alphar.tau0delta1(tau, delta) + self.alphar.tau0delta0(tau, delta) + delta.ln()
    }
    fn Jdelta(&self, tau: f64, delta: f64) -> f64 {
        1.0 + 2.0 * self.alphar.tau0delta1(tau, delta) + self.alphar.tau0delta2(tau, delta)
    }
    fn Kdelta(&self, tau: f64, delta: f64) -> f64 {
        (2.0 * self.alphar.tau0delta1(tau, delta) + self.alphar.tau0delta2(tau, delta) + 1.0)
            / delta
    }
}
#[pymethods]
#[allow(non_snake_case)]
impl Helmholtz {
    /// read helmholtz equation of state from fluid.json file
    #[new]
    pub fn read_json(path: &str) -> anyhow::Result<Helmholtz> {
        let mut file_json: Option<std::fs::File>;
        file_json = match std::fs::File::open(std::path::Path::new(path)) {
            Ok(file) => Some(file),
            Err(_) => None,
        };
        if file_json.is_none() {
            file_json = match std::fs::File::open(
                std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("res")
                    .join(path),
            ) {
                Ok(file) => Some(file),
                Err(_) => None,
            };
        };
        let mut file_json = file_json.ok_or(anyhow!(HelmholtzErr::NoJson))?;
        let mut str_json = String::new();
        let _ = std::io::Read::read_to_string(&mut file_json, &mut str_json);
        let eos: Option<Helmholtz> = match serde_json::from_str(&str_json) {
            Ok(eos) => Some(eos),
            Err(_) => None,
        };
        let eos = eos.ok_or(anyhow!(HelmholtzErr::NoHelmholtz))?;
        Ok(eos)
    }
    pub fn td_unchecked(&mut self, T: f64, rho: f64) {
        self.phase = Phase::One { T, rho }
    }
    pub fn t_flash(&mut self, Ts: f64) -> anyhow::Result<()> {
        let tau = self.Tc / Ts;
        let delta_v0 = self.rhovs.calc(Ts, self.Tc, 1.0);
        let delta_l0 = self.rhols.calc(Ts, self.Tc, 1.0);
        // Reference [AKASAKA_2008]
        // Newton-Raphson algorithm
        let (mut Jv, mut DJv, mut Kv, mut DKv): (f64, f64, f64, f64);
        let (mut Jl, mut DJl, mut Kl, mut DKl): (f64, f64, f64, f64);
        let mut delta_v = 1.0;
        let mut delta_l = 1.0;
        let mut Delta: f64;
        let mut is_conv: bool = false;
        'outer: for factor in 0..10 {
            delta_v = delta_v0 * (1.0 - (factor as f64) / 10.0);
            delta_l = delta_l0 * (1.0 + (factor as f64) / 10.0);
            'inner: loop {
                Jv = self.J(tau, delta_v);
                Jl = self.J(tau, delta_l);
                Kv = self.K(tau, delta_v);
                Kl = self.K(tau, delta_l);
                // The following convergence condition was used:
                if ((Kv - Kl).abs() + (Jv - Jl).abs()) < 1E-8 {
                    is_conv = true;
                    break 'outer; // convergence
                }
                DJv = self.Jdelta(tau, delta_v);
                DJl = self.Jdelta(tau, delta_l);
                // Jdelta = DP_DD_T /RT
                if DJv.is_sign_negative() || DJl.is_sign_negative() {
                    break 'inner; //divergence
                }
                DKv = self.Kdelta(tau, delta_v);
                DKl = self.Kdelta(tau, delta_l);
                Delta = DJv * DKl - DJl * DKv;
                delta_v += ((Kv - Kl) * DJl - (Jv - Jl) * DKl) / Delta;
                delta_l += ((Kv - Kl) * DJv - (Jv - Jl) * DKv) / Delta;
                // make sure rhog != NAN and rhol != NAN
                if delta_v.is_nan() || delta_l.is_nan() {
                    break 'inner; //divergence
                }
            }
        }
        if is_conv {
            self.phase = Phase::Two {
                Ts,
                rhov: delta_v * self.rhoc,
                rhol: delta_l * self.rhoc,
                x: 1.0,
            };
            Ok(())
        } else {
            Err(anyhow!(HelmholtzErr::NotConvForT))
        }
    }
    pub fn td_flash(&mut self, T: f64, rho: f64) -> anyhow::Result<()> {
        if T >= self.Tc
            || rho <= 0.85 * self.rhovs.calc(T, self.Tc, self.rhoc)
            || rho >= 1.05 * self.rhols.calc(T, self.Tc, self.rhoc)
        {
            self.phase = Phase::One { T, rho };
            Ok(()) // one phase
        } else if let Err(why) = self.t_flash(T) {
            Err(why) // Self::t_flash diverge
        } else {
            match self.phase {
                Phase::Two { rhov, rhol, .. } => {
                    if rho < rhov || rho > rhol {
                        self.phase = Phase::One { T, rho };
                        Ok(()) // one phase
                    } else {
                        self.phase = Phase::Two {
                            Ts: T,
                            rhov,
                            rhol,
                            x: (1.0 / rho - 1.0 / rhol) / (1.0 / rhov - 1.0 / rhol),
                        };
                        Ok(()) // two phase
                    }
                }
                _ => Err(anyhow!(HelmholtzErr::NotConvForTD)), // shit-code in td_flash
            }
        }
    }
    pub fn tp_flash(&mut self, T: f64, p: f64) -> anyhow::Result<()> {
        let phase: char;
        let ps = self.ps.calc(T, self.Tc, self.pc);
        if T >= self.Tc {
            phase = 's'; // supercritical region
        } else if p < 0.95 * ps {
            phase = 'v' // vapor region
        } else if p > 1.05 * ps {
            phase = 'l'; // liquid region
        } else if let Err(why) = self.t_flash(T) {
            return Err(why);
        } else {
            let ps = self.p_s().unwrap();
            if p < ps {
                phase = 'v';
            } else {
                phase = 'l';
            }
        }
        let mut rho = rho_pr(T, p, self.Tc, self.pc, self.R, self.omega);
        if phase == 'v' {
            rho = rho.min(self.rhovs.calc(T, self.Tc, self.rhoc));
        }
        if phase == 'l' {
            rho = rho.max(self.rhols.calc(T, self.Tc, self.rhoc));
        }
        let mut p_diff: f64;
        loop {
            p_diff = self.calc_p(T, rho) - p;
            if (p_diff / p).abs() < 1E-9 {
                break;
            }
            rho -= p_diff / self.calc_Dp_Drho_T(T, rho);
        }
        self.phase = Phase::One { T, rho };
        Ok(())
    }
    pub fn T(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { T, .. } => Ok(T),
            Phase::Two { .. } => Err(anyhow!(HelmholtzErr::NotInTwoPhase)),
        }
    }
    pub fn p(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { T, rho } => Ok(self.calc_p(T, rho)),
            Phase::Two { .. } => Err(anyhow!(HelmholtzErr::NotInTwoPhase)),
        }
    }
    pub fn rho(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { rho, .. } => Ok(rho),
            Phase::Two { .. } => Err(anyhow!(HelmholtzErr::NotInTwoPhase)),
        }
    }
    pub fn cv(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { T, rho } => Ok(self.calc_cv(T, rho)),
            Phase::Two { .. } => Err(anyhow!(HelmholtzErr::NotInTwoPhase)),
        }
    }
    pub fn cp(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { T, rho } => Ok(self.calc_cp(T, rho)),
            Phase::Two { .. } => Err(anyhow!(HelmholtzErr::NotInTwoPhase)),
        }
    }
    pub fn w(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { T, rho } => Ok(self.calc_w(T, rho)),
            Phase::Two { .. } => Err(anyhow!(HelmholtzErr::NotInTwoPhase)),
        }
    }
    pub fn h(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { T, rho } => Ok(self.calc_h(T, rho)),
            Phase::Two { Ts, rhov, rhol, x } => {
                Ok(x * self.calc_h(Ts, rhov) + (1.0 - x) * self.calc_h(Ts, rhol))
            }
        }
    }
    pub fn s(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { T, rho } => Ok(self.calc_s(T, rho)),
            Phase::Two { Ts, rhov, rhol, x } => {
                Ok(x * self.calc_s(Ts, rhov) + (1.0 - x) * self.calc_s(Ts, rhol))
            }
        }
    }
    pub fn T_s(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { .. } => Err(anyhow!(HelmholtzErr::NotInOnePhase)),
            Phase::Two { Ts, .. } => Ok(Ts),
        }
    }
    pub fn p_s(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { .. } => Err(anyhow!(HelmholtzErr::NotInOnePhase)),
            Phase::Two { Ts, rhov, rhol, x } => {
                Ok(x * self.calc_p(Ts, rhov) + (1.0 - x) * self.calc_p(Ts, rhol))
            }
        }
    }
    pub fn rho_v(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { .. } => Err(anyhow!(HelmholtzErr::NotInOnePhase)),
            Phase::Two { rhov, .. } => Ok(rhov),
        }
    }
    pub fn rho_l(&self) -> anyhow::Result<f64> {
        match self.phase {
            Phase::One { .. } => Err(anyhow!(HelmholtzErr::NotInOnePhase)),
            Phase::Two { rhol, .. } => Ok(rhol),
        }
    }
}
