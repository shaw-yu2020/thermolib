use super::ancillary_equations::SaturatedLiquidDensityEquation;
use super::ancillary_equations::SaturatedVaporDensityEquation;
use super::ancillary_equations::VaporPressureEquation;
use super::ideal_helmholtz_equation::IdealHelmholtzEquation;
use super::residual_helmholtz_equation::ResidualHelmholtzEquation;
use super::AlphaDD;
use super::ThermoProp;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Read;
use std::path::Path;
///
/// HelmholtzPure 状态方程模型
///
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug)]
pub struct HelmholtzPure {
    Tc: f64,    // 临界温度
    Pc: f64,    // 临界压力
    Dc: f64,    // 临界密度
    R: f64,     // 气体常数
    M: f64,     // 摩尔质量
    omega: f64, // 偏心因子
    alpha0: IdealHelmholtzEquation,
    alphar: ResidualHelmholtzEquation,
    ps: VaporPressureEquation,
    rhogs: SaturatedVaporDensityEquation,
    rhols: SaturatedLiquidDensityEquation,
}
pub fn read_json(path_json: &str) -> HelmholtzPure {
    let path = Path::new(path_json);
    let mut file = match File::open(&path) {
        Ok(file) => file,
        Err(why) => {
            print!("couldn't open {}: {:?}\n", path.display(), why);
            let path_thermolib = Path::new(env!("CARGO_MANIFEST_DIR"));
            let path = path_thermolib.join("res").join(path_json);
            print!("search {}\n", path.display());
            match File::open(&path) {
                Ok(file) => file,
                Err(why) => panic!("couldn't open {}: {:?}", path.display(), why),
            }
        }
    };
    let mut str_json = String::new();
    match file.read_to_string(&mut str_json) {
        Ok(_) => print!("{} contains:\n{}\n", path.display(), str_json),
        Err(why) => panic!("couldn't read {}: {:?}", path.display(), why),
    }
    serde_json::from_str(&str_json).unwrap()
}
impl HelmholtzPure {
    #[allow(non_snake_case)]
    fn calc(&self, tp: ThermoProp, T: f64, D: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = D / self.Dc;
        match tp {
            ThermoProp::T => T,
            ThermoProp::D => D,
            ThermoProp::Z => 1.0 + self.alphar.calc(AlphaDD::D01, tau, delta),
            ThermoProp::P => (D * self.R * T) * (1.0 + self.alphar.calc(AlphaDD::D01, tau, delta)),
            ThermoProp::CV => {
                self.R
                    * (-tau.powi(2) * self.alpha0.calc(AlphaDD::D20, tau, delta, self.Tc)
                        - self.alphar.calc(AlphaDD::D20, tau, delta))
            }
            ThermoProp::CP => {
                self.R
                    * ((-tau.powi(2) * self.alpha0.calc(AlphaDD::D20, tau, delta, self.Tc)
                        - self.alphar.calc(AlphaDD::D20, tau, delta))
                        + (1.0 + self.alphar.calc(AlphaDD::D01, tau, delta)
                            - self.alphar.calc(AlphaDD::D11, tau, delta))
                        .powi(2)
                            / (1.0
                                + 2.0 * self.alphar.calc(AlphaDD::D01, tau, delta)
                                + self.alphar.calc(AlphaDD::D02, tau, delta)))
            }
            ThermoProp::W => (if self.R < 10.0 {
                self.R / self.M
            } else {
                self.R
            } * T
                * ((1.0
                    + 2.0 * self.alphar.calc(AlphaDD::D01, tau, delta)
                    + self.alphar.calc(AlphaDD::D02, tau, delta))
                    + (1.0 + self.alphar.calc(AlphaDD::D01, tau, delta)
                        - self.alphar.calc(AlphaDD::D11, tau, delta))
                    .powi(2)
                        / (-tau.powi(2) * self.alpha0.calc(AlphaDD::D20, tau, delta, self.Tc)
                            - self.alphar.calc(AlphaDD::D20, tau, delta))))
            .sqrt(),
            ThermoProp::S => {
                self.R
                    * (tau * self.alpha0.calc(AlphaDD::D10, tau, delta, self.Tc)
                        + self.alphar.calc(AlphaDD::D10, tau, delta)
                        - self.alpha0.calc(AlphaDD::D00, tau, delta, self.Tc)
                        - self.alphar.calc(AlphaDD::D00, tau, delta))
            }
            ThermoProp::U => {
                (self.R * T)
                    * (tau * self.alpha0.calc(AlphaDD::D10, tau, delta, self.Tc)
                        + self.alphar.calc(AlphaDD::D10, tau, delta))
            }
            ThermoProp::H => {
                (self.R * T)
                    * (1.0
                        + tau * self.alpha0.calc(AlphaDD::D10, tau, delta, self.Tc)
                        + self.alphar.calc(AlphaDD::D10, tau, delta)
                        + self.alphar.calc(AlphaDD::D01, tau, delta))
            }
            ThermoProp::A => {
                (self.R * T)
                    * (self.alpha0.calc(AlphaDD::D00, tau, delta, self.Tc)
                        + self.alphar.calc(AlphaDD::D00, tau, delta))
            }
            ThermoProp::G => {
                (self.R * T)
                    * (1.0
                        + self.alpha0.calc(AlphaDD::D00, tau, delta, self.Tc)
                        + self.alphar.calc(AlphaDD::D00, tau, delta)
                        + self.alphar.calc(AlphaDD::D01, tau, delta))
            }
        }
    }
}
