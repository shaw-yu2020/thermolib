use super::{LiquidMetal, LiquidMetalErr, Metals};
use anyhow::anyhow;
use lazy_static::lazy_static;
use std::collections::HashMap;
#[allow(non_snake_case)]
impl LiquidMetal {
    /// calculate viscosity of liquid metals at 0.1 MPa
    pub fn eta(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(params) = MAP.get(&self.metal) {
            if T < params.Tmin {
                Err(anyhow!(LiquidMetalErr::TisTooMin))
            } else if T > params.Tmax {
                Err(anyhow!(LiquidMetalErr::TisTooMax))
            } else {
                Ok(10.0_f64.powf(params.a0 + params.a1 / T) * params.eta0)
            }
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
}
#[allow(non_snake_case)]
struct EtaParams {
    Tmin: f64,
    Tmax: f64,
    a0: f64,
    a1: f64,
    eta0: f64,
}
lazy_static! {
    static ref MAP: HashMap<Metals, EtaParams> = HashMap::from([]);
}
