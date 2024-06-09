mod alpha_i;
use alpha_i::IdealHelmholtz;
mod alpha_r;
use alpha_r::ResidualHelmholtz;
mod anc_eqn;
use anc_eqn::{PsEqn, RholEqn, RhovEqn};
mod rho_ini;
use rho_ini::rho_pr;
/// Helmholtz EOS
mod alpha;
pub use alpha::Helmholtz;
/// HelmholtzError
use thiserror::Error;
#[derive(Debug, Error)]
pub enum HelmholtzErr {
    #[error("t_flash diverge")]
    NotConvForT,
    #[error("td_flash diverge")]
    NotConvForTD,
    #[error("not in one phase")]
    NotInOnePhase,
    #[error("not in two phase")]
    NotInTwoPhase,
    #[error("no fluid.json")]
    NoJson,
    #[error("no helmholtz")]
    NoHelmholtz,
}
#[allow(non_snake_case)]
pub enum Phase {
    One {
        T: f64,
        rho: f64,
    },
    Two {
        Ts: f64,
        rhov: f64,
        rhol: f64,
        x: f64,
    },
}
impl Default for Phase {
    fn default() -> Self {
        Phase::One { T: 1.0, rho: 0.0 }
    }
}
