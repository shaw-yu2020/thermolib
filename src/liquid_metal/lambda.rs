use super::{LiquidMetal, LiquidMetalErr, Metals};
use anyhow::anyhow;
use lazy_static::lazy_static;
use std::collections::HashMap;
#[allow(non_snake_case)]
impl LiquidMetal {
    /// calculate thermal conductivity of liquid metals at 0.1 MPa
    pub fn lambda(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(params) = MAP.get(&self.metal) {
            if T < params.Tmin {
                Err(anyhow!(LiquidMetalErr::TisTooMin))
            } else if T > params.Tmax {
                Err(anyhow!(LiquidMetalErr::TisTooMax))
            } else {
                Ok(params.d0 + params.d1 * (T - params.Tm) + params.d2 * (T - params.Tm).powi(2))
            }
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
}
#[allow(non_snake_case)]
struct LambdaParams {
    Tm: f64,
    Tmin: f64,
    Tmax: f64,
    d0: f64,
    d1: f64,
    d2: f64,
}
lazy_static! {
    static ref MAP: HashMap<Metals, LambdaParams> = HashMap::from([
        (
            Metals::Ti,
            LambdaParams {
                Tm: 1941.0,
                Tmin: 1941.0,
                Tmax: 5000.0,
                d0: 30.693,
                d1: 12.294E-3,
                d2: -11.982E-7
            },
        ),
        (
            Metals::V,
            LambdaParams {
                Tm: 2183.0,
                Tmin: 2183.0,
                Tmax: 3900.0,
                d0: 42.045,
                d1: 14.007E-3,
                d2: -26.085E-7
            },
        ),
        (
            Metals::Cr,
            LambdaParams {
                Tm: 2180.0,
                Tmin: 2186.0,
                Tmax: 2198.0,
                d0: 0.0,
                d1: 0.0,
                d2: 0.0
            },
        ),
        (
            Metals::Zr,
            LambdaParams {
                Tm: 2128.0,
                Tmin: 2128.0,
                Tmax: 4275.0,
                d0: 38.151,
                d1: 15.074E-3,
                d2: -13.172E-7
            },
        ),
        (
            Metals::Nb,
            LambdaParams {
                Tm: 2742.0,
                Tmin: 2742.0,
                Tmax: 4450.0,
                d0: 62.68,
                d1: 19.386E-3,
                d2: -30.596E-7
            },
        ),
        (
            Metals::Mo,
            LambdaParams {
                Tm: 2896.0,
                Tmin: 2896.0,
                Tmax: 4500.0,
                d0: 71.832,
                d1: 23.872E-3,
                d2: -34.946E-7
            },
        ),
        (
            Metals::Hf,
            LambdaParams {
                Tm: 2500.0,
                Tmin: 2500.0,
                Tmax: 3500.0,
                d0: 37.891,
                d1: 14.08E-3,
                d2: 0.0
            },
        ),
        (
            Metals::Ta,
            LambdaParams {
                Tm: 3293.0,
                Tmin: 3293.0,
                Tmax: 6900.0,
                d0: 62.201,
                d1: 16.493E-3,
                d2: -13.04E-7
            },
        ),
        (
            Metals::W,
            LambdaParams {
                Tm: 3695.0,
                Tmin: 3695.0,
                Tmax: 5800.0,
                d0: 65.26,
                d1: 18.595E-3,
                d2: -20.217E-7
            },
        )
    ]);
}
