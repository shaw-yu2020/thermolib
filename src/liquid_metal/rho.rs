use super::{LiquidMetalErr, Metals};
use anyhow::anyhow;
use std::collections::HashMap;
use std::sync::OnceLock;
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
pub fn metals2rhoparams() -> &'static HashMap<Metals, RhoParams> {
    static METALS_TO_RHOPARAMS: OnceLock<HashMap<Metals, RhoParams>> = OnceLock::new();
    METALS_TO_RHOPARAMS.get_or_init(|| {
        let mut hm = HashMap::new();
        hm.insert(
            Metals::Al,
            RhoParams {
                Tm: 933.47,
                Tmin: 933.0,
                Tmax: 1190.0,
                c0: 2377.23,
                c1: -0.311,
            },
        );
        hm.insert(
            Metals::Si,
            RhoParams {
                Tm: 1687.0,
                Tmin: 1687.0,
                Tmax: 2000.0,
                c0: 2550.0,
                c1: -0.264,
            },
        );
        hm.insert(
            Metals::Ti,
            RhoParams {
                Tm: 1941.0,
                Tmin: 1941.0,
                Tmax: 3520.0,
                c0: 4222.1,
                c1: -0.3952,
            },
        );
        hm.insert(
            Metals::V,
            RhoParams {
                Tm: 2183.0,
                Tmin: 2183.0,
                Tmax: 4500.0,
                c0: 5517.0,
                c1: -0.5895,
            },
        );
        hm.insert(
            Metals::Cr,
            RhoParams {
                Tm: 2180.0,
                Tmin: 2186.0,
                Tmax: 2503.0,
                c0: 6097.1,
                c1: -0.6536,
            },
        );
        hm.insert(
            Metals::Fe,
            RhoParams {
                Tm: 1811.0,
                Tmin: 1809.0,
                Tmax: 2480.0,
                c0: 7034.96,
                c1: -0.926,
            },
        );
        hm.insert(
            Metals::Co,
            RhoParams {
                Tm: 1768.0,
                Tmin: 1768.0,
                Tmax: 2500.0,
                c0: 7827.0,
                c1: -0.936,
            },
        );
        hm.insert(
            Metals::Ni,
            RhoParams {
                Tm: 1728.0,
                Tmin: 1728.0,
                Tmax: 2500.0,
                c0: 7861.0,
                c1: -0.988,
            },
        );
        hm.insert(
            Metals::Cu,
            RhoParams {
                Tm: 1357.77,
                Tmin: 1356.0,
                Tmax: 2500.0,
                c0: 7997.0,
                c1: -0.819,
            },
        );
        hm.insert(
            Metals::Zn,
            RhoParams {
                Tm: 692.677,
                Tmin: 692.0,
                Tmax: 910.0,
                c0: 6559.0,
                c1: -0.884,
            },
        );
        hm.insert(
            Metals::Ga,
            RhoParams {
                Tm: 302.914,
                Tmin: 303.0,
                Tmax: 1500.0,
                c0: 6077.0,
                c1: -0.611,
            },
        );
        hm.insert(
            Metals::Zr,
            RhoParams {
                Tm: 2128.0,
                Tmin: 2128.0,
                Tmax: 4100.0,
                c0: 6100.0,
                c1: -0.242,
            },
        );
        hm.insert(
            Metals::Nb,
            RhoParams {
                Tm: 2742.0,
                Tmin: 2742.0,
                Tmax: 5848.0,
                c0: 7664.0,
                c1: -0.2943,
            },
        );
        hm.insert(
            Metals::Mo,
            RhoParams {
                Tm: 2896.0,
                Tmin: 2896.0,
                Tmax: 5914.0,
                c0: 9062.6,
                c1: -0.3947,
            },
        );
        hm.insert(
            Metals::Ag,
            RhoParams {
                Tm: 1234.93,
                Tmin: 1235.0,
                Tmax: 1600.0,
                c0: 9294.0,
                c1: -0.877,
            },
        );
        hm.insert(
            Metals::Cd,
            RhoParams {
                Tm: 594.219,
                Tmin: 594.0,
                Tmax: 833.0,
                c0: 8008.0,
                c1: -1.251,
            },
        );
        hm.insert(
            Metals::In,
            RhoParams {
                Tm: 429.748,
                Tmin: 430.0,
                Tmax: 1100.0,
                c0: 7022.0,
                c1: -0.762,
            },
        );
        hm.insert(
            Metals::Sn,
            RhoParams {
                Tm: 505.08,
                Tmin: 506.0,
                Tmax: 1950.0,
                c0: 6979.0,
                c1: -0.652,
            },
        );
        hm.insert(
            Metals::Sb,
            RhoParams {
                Tm: 899.0,
                Tmin: 900.0,
                Tmax: 1300.0,
                c0: 6467.0,
                c1: -0.608,
            },
        );
        hm.insert(
            Metals::Hf,
            RhoParams {
                Tm: 2500.0,
                Tmin: 2500.0,
                Tmax: 4981.0,
                c0: 11902.6,
                c1: -0.6704,
            },
        );
        hm.insert(
            Metals::Ta,
            RhoParams {
                Tm: 3293.0,
                Tmin: 3293.0,
                Tmax: 6400.0,
                c0: 14977.5,
                c1: -0.6802,
            },
        );
        hm.insert(
            Metals::W,
            RhoParams {
                Tm: 3695.0,
                Tmin: 3695.0,
                Tmax: 5818.0,
                c0: 17146.4,
                c1: -0.6769,
            },
        );
        hm.insert(
            Metals::Tl,
            RhoParams {
                Tm: 576.7,
                Tmin: 576.0,
                Tmax: 1200.0,
                c0: 11233.0,
                c1: -1.2,
            },
        );
        hm.insert(
            Metals::Pb,
            RhoParams {
                Tm: 600.61,
                Tmin: 601.0,
                Tmax: 2000.0,
                c0: 10656.0,
                c1: -1.239,
            },
        );
        hm.insert(
            Metals::Bi,
            RhoParams {
                Tm: 544.55,
                Tmin: 545.0,
                Tmax: 1500.0,
                c0: 10028.0,
                c1: -1.213,
            },
        );
        hm
    })
}
