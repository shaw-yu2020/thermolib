use super::{LiquidMetalErr, Metals};
use anyhow::anyhow;
use lazy_static::lazy_static;
use std::collections::HashMap;
#[allow(non_snake_case)]
pub struct LambdaParams {
    Tm: f64,
    Tmin: f64,
    Tmax: f64,
    d0: f64,
    d1: f64,
    d2: f64,
}
#[allow(non_snake_case)]
impl LambdaParams {
    pub fn calc(&self, T: f64) -> anyhow::Result<f64> {
        if T < self.Tmin {
            Err(anyhow!(LiquidMetalErr::TisTooMin))
        } else if T > self.Tmax {
            Err(anyhow!(LiquidMetalErr::TisTooMax))
        } else {
            Ok(self.d0 + self.d1 * (T - self.Tm) + self.d2 * (T - self.Tm).powi(2))
        }
    }
}
lazy_static! {
    pub static ref METALS_TO_LAMBDAPARAMS: HashMap<Metals, LambdaParams> = HashMap::from([
        (
            Metals::Si,
            LambdaParams {
                Tm: 1687.0,
                Tmin: 1690.0,
                Tmax: 1945.0,
                d0: 54.70218,
                d1: 1.53E-3,
                d2: 0.0
            }
        ),
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
            Metals::Fe,
            LambdaParams {
                Tm: 1811.0,
                Tmin: 1815.0,
                Tmax: 2050.0,
                d0: 36.349,
                d1: 9.6207E-3,
                d2: 0.0
            }
        ),
        (
            Metals::Co,
            LambdaParams {
                Tm: 1768.15,
                Tmin: 1769.0,
                Tmax: 1903.0,
                d0: 29.49359,
                d1: 87.81E-3,
                d2: 0.0
            }
        ),
        (
            Metals::Ni,
            LambdaParams {
                Tm: 1728.0,
                Tmin: 1730.0,
                Tmax: 2000.0,
                d0: 54.182,
                d1: 20.97E-3,
                d2: 0.0
            }
        ),
        (
            Metals::Cu,
            LambdaParams {
                Tm: 1357.77,
                Tmin: 1358.0,
                Tmax: 1700.0,
                d0: 150.49,
                d1: 70.41E-3,
                d2: 0.0
            }
        ),
        (
            Metals::Ga,
            LambdaParams {
                Tm: 302.914,
                Tmin: 303.0,
                Tmax: 850.0,
                d0: 28.403,
                d1: 71.896E-3,
                d2: 0.0
            }
        ),
        (
            Metals::Ge,
            LambdaParams {
                Tm: 1210.4,
                Tmin: 1212.0,
                Tmax: 1473.0,
                d0: 45.55252,
                d1: 24.09E-3,
                d2: 0.0
            }
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
            Metals::In,
            LambdaParams {
                Tm: 429.748,
                Tmin: 430.0,
                Tmax: 1300.0,
                d0: 36.493,
                d1: 29.185E-3,
                d2: 0.0
            }
        ),
        (
            Metals::Sn,
            LambdaParams {
                Tm: 505.8,
                Tmin: 507.0,
                Tmax: 2000.0,
                d0: 28.037,
                d1: 23.397E-3,
                d2: 0.0
            }
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
        ),
        (
            Metals::Pb,
            LambdaParams {
                Tm: 600.61,
                Tmin: 602.0,
                Tmax: 1150.0,
                d0: 16.093,
                d1: 7.8526E-3,
                d2: 0.0
            }
        ),
        (
            Metals::Bi,
            LambdaParams {
                Tm: 544.55,
                Tmin: 545.0,
                Tmax: 1110.0,
                d0: 13.19939,
                d1: 11.47E-3,
                d2: 0.0
            }
        )
    ]);
}
