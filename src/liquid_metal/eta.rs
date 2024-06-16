use super::{LiquidMetalErr, Metals};
use anyhow::anyhow;
use lazy_static::lazy_static;
use std::collections::HashMap;
#[allow(non_snake_case)]
pub struct EtaParams {
    Tmin: f64,
    Tmax: f64,
    a0: f64,
    a1: f64,
    eta0: f64,
}
#[allow(non_snake_case)]
impl EtaParams {
    pub fn calc(&self, T: f64) -> anyhow::Result<f64> {
        if T < self.Tmin {
            Err(anyhow!(LiquidMetalErr::TisTooMin))
        } else if T > self.Tmax {
            Err(anyhow!(LiquidMetalErr::TisTooMax))
        } else {
            Ok(10.0_f64.powf(self.a0 + self.a1 / T) * self.eta0)
        }
    }
}
lazy_static! {
    pub static ref METALS_TO_ETAPARAMS: HashMap<Metals, EtaParams> = HashMap::from([
        (
            Metals::Al,
            EtaParams {
                Tmin: 933.0,
                Tmax: 1270.0,
                a0: -0.7324,
                a1: 803.49,
                eta0: 1.0,
            },
        ),
        (
            Metals::Fe,
            EtaParams {
                Tmin: 1809.0,
                Tmax: 2480.0,
                a0: -0.7209,
                a1: 2694.85,
                eta0: 1.0
            }
        ),
        (
            Metals::Cu,
            EtaParams {
                Tmin: 1356.0,
                Tmax: 1970.0,
                a0: -0.422,
                a1: 1393.4,
                eta0: 1.0
            }
        ),
        (
            Metals::Sn,
            EtaParams {
                Tmin: 506.0,
                Tmax: 1280.0,
                a0: -0.408,
                a1: 343.4,
                eta0: 1.0
            }
        )
    ]);
}
