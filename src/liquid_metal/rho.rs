use super::{LiquidMetal, LiquidMetalErr, Metals};
use anyhow::anyhow;
use lazy_static::lazy_static;
use std::collections::HashMap;
#[allow(non_snake_case)]
impl LiquidMetal {
    /// calculate density of liquid metals at 0.1 MPa
    pub fn rho(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(params) = MAP.get(&self.metal) {
            if T < params.Tmin {
                Err(anyhow!(LiquidMetalErr::TisTooMin))
            } else if T > params.Tmax {
                Err(anyhow!(LiquidMetalErr::TisTooMax))
            } else {
                Ok(params.c0 + params.c1 * (T - params.Tm))
            }
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
}
#[allow(non_snake_case)]
struct RhoParams {
    Tm: f64,
    Tmin: f64,
    Tmax: f64,
    c0: f64,
    c1: f64,
}
lazy_static! {
    static ref MAP: HashMap<Metals, RhoParams> = HashMap::from([
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
