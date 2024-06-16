use super::{LiquidMetalErr, Metals};
use anyhow::anyhow;
use lazy_static::lazy_static;
use std::collections::HashMap;
#[allow(non_snake_case)]
pub struct RhoParams {
    Tm: f64,
    Tmin: f64,
    Tmax: f64,
    c0: f64,
    c1: f64,
}
#[allow(non_snake_case)]
impl RhoParams {
    pub fn calc(&self, T: f64) -> anyhow::Result<f64> {
        if T < self.Tmin {
            Err(anyhow!(LiquidMetalErr::TisTooMin))
        } else if T > self.Tmax {
            Err(anyhow!(LiquidMetalErr::TisTooMax))
        } else {
            Ok(self.c0 + self.c1 * (T - self.Tm))
        }
    }
}
lazy_static! {
    pub static ref METALS_TO_RHOPARAMS: HashMap<Metals, RhoParams> = HashMap::from([
        (
            Metals::Al,
            RhoParams {
                Tm: 933.47,
                Tmin: 933.0,
                Tmax: 1190.0,
                c0: 2377.23,
                c1: -0.311
            }
        ),
        (
            Metals::Ti,
            RhoParams {
                Tm: 1941.0,
                Tmin: 1941.0,
                Tmax: 3520.0,
                c0: 4222.1,
                c1: -0.3952,
            },
        ),
        (
            Metals::V,
            RhoParams {
                Tm: 2183.0,
                Tmin: 2183.0,
                Tmax: 4500.0,
                c0: 5517.0,
                c1: -0.5895
            }
        ),
        (
            Metals::Cr,
            RhoParams {
                Tm: 2180.0,
                Tmin: 2186.0,
                Tmax: 2503.0,
                c0: 6097.1,
                c1: -0.6536
            }
        ),
        (
            Metals::Fe,
            RhoParams {
                Tm: 1811.0,
                Tmin: 1809.0,
                Tmax: 2480.0,
                c0: 7034.96,
                c1: -0.926
            }
        ),
        (
            Metals::Zr,
            RhoParams {
                Tm: 2128.0,
                Tmin: 2128.0,
                Tmax: 4100.0,
                c0: 6100.0,
                c1: -0.242
            }
        ),
        (
            Metals::Nb,
            RhoParams {
                Tm: 2742.0,
                Tmin: 2742.0,
                Tmax: 5848.0,
                c0: 7664.0,
                c1: -0.2943
            }
        ),
        (
            Metals::Mo,
            RhoParams {
                Tm: 2896.0,
                Tmin: 2896.0,
                Tmax: 5914.0,
                c0: 9062.6,
                c1: -0.3947
            }
        ),
        (
            Metals::Hf,
            RhoParams {
                Tm: 2500.0,
                Tmin: 2500.0,
                Tmax: 4981.0,
                c0: 11902.6,
                c1: -0.6704
            }
        ),
        (
            Metals::Ta,
            RhoParams {
                Tm: 3293.0,
                Tmin: 3293.0,
                Tmax: 6400.0,
                c0: 14977.5,
                c1: -0.6802
            }
        ),
        (
            Metals::W,
            RhoParams {
                Tm: 3695.0,
                Tmin: 3695.0,
                Tmax: 5818.0,
                c0: 17146.4,
                c1: -0.6769
            }
        )
    ]);
}
