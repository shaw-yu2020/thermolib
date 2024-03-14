use super::ancillary_equations::SaturatedLiquidDensityEquation;
use super::ancillary_equations::SaturatedVaporDensityEquation;
use super::ancillary_equations::VaporPressureEquation;
use super::real_helmholtz_equation::RealHelmholtzEquation;
use super::ThermoProp;
use crate::pubtraits::Flash;
use crate::pubtraits::MyErr;
use crate::pubtraits::Prop;
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
    alpha: RealHelmholtzEquation,
    ps: VaporPressureEquation,
    rhogs: SaturatedVaporDensityEquation,
    rhols: SaturatedLiquidDensityEquation,
    omega: f64, // 偏心因子
    #[serde(skip, default = "default_phase")]
    phase: Phase, // 记录相态
    #[serde(skip)]
    T: f64, // 记录温度
}
#[derive(Debug)]
enum Phase {
    SINGLE { rho: f64 },                     // 密度
    SATVAP { rhog: f64 },                    // 饱和气相密度
    SATLIQ { rhol: f64 },                    // 饱和液相密度
    DOUBLE { rhog: f64, rhol: f64, x: f64 }, // 饱和气相密度 饱和液相密度 干度
}
fn default_phase() -> Phase {
    Phase::SINGLE { rho: 1.0 }
}
pub fn read_json(path: &str) -> Result<HelmholtzPure, MyErr> {
    let mut str_json = String::new();
    let _ = match File::open(&Path::new(path)) {
        Ok(file) => file,
        Err(_) => match File::open(&Path::new(env!("CARGO_MANIFEST_DIR")).join("res").join(path)) {
            Ok(file) => file,
            Err(_) => return Err(MyErr::new(&format!("couldn't find {}", path))),
        },
    }
    .read_to_string(&mut str_json);
    match serde_json::from_str(&str_json) {
        Ok(hp) => Ok(hp),
        Err(_) => Err(MyErr::new(&format!("no alpha(HelmholtzPure) in {}", path))),
    }
}
impl Flash for HelmholtzPure {
    #[allow(non_snake_case)]
    fn td_flash(&mut self, T: f64, rho: f64) -> Result<(), MyErr> {
        self.T = T;
        self.phase = Phase::SINGLE { rho };
        Ok(())
    }
}
impl Prop for HelmholtzPure {
    fn T(&self) -> Result<f64, MyErr> {
        Ok(self.T)
    }
    fn rho(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(rho),
            Phase::SATVAP { rhog } => Ok(rhog),
            Phase::SATLIQ { rhol } => Ok(rhol),
            Phase::DOUBLE { x, rhog, rhol } => Ok(1.0 / (x / rhog + (1.0 - x) / rhol)),
        }
    }
    fn p(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.alpha.calc(ThermoProp::P, self.T, rho)),
            Phase::SATVAP { rhog } => Ok(self.alpha.calc(ThermoProp::P, self.T, rhog)),
            Phase::SATLIQ { rhol } => Ok(self.alpha.calc(ThermoProp::P, self.T, rhol)),
            Phase::DOUBLE { x, rhog, rhol } => Ok(x * self.alpha.calc(ThermoProp::P, self.T, rhog)
                + (1.0 - x) * self.alpha.calc(ThermoProp::P, self.T, rhol)),
        }
    }
    fn cv(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.alpha.calc(ThermoProp::CV, self.T, rho)),
            Phase::SATVAP { rhog } => Ok(self.alpha.calc(ThermoProp::CV, self.T, rhog)),
            Phase::SATLIQ { rhol } => Ok(self.alpha.calc(ThermoProp::CV, self.T, rhol)),
            Phase::DOUBLE { x, rhog, rhol } => {
                Ok(x * self.alpha.calc(ThermoProp::CV, self.T, rhog)
                    + (1.0 - x) * self.alpha.calc(ThermoProp::CV, self.T, rhol))
            }
        }
    }
    fn cp(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.alpha.calc(ThermoProp::CP, self.T, rho)),
            Phase::SATVAP { rhog } => Ok(self.alpha.calc(ThermoProp::CP, self.T, rhog)),
            Phase::SATLIQ { rhol } => Ok(self.alpha.calc(ThermoProp::CP, self.T, rhol)),
            Phase::DOUBLE { x, rhog, rhol } => {
                Ok(x * self.alpha.calc(ThermoProp::CP, self.T, rhog)
                    + (1.0 - x) * self.alpha.calc(ThermoProp::P, self.T, rhol))
            }
        }
    }
    fn w(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.alpha.calc(ThermoProp::W, self.T, rho)),
            _ => Err(MyErr::new(&format!("no speed of sound in double phase"))),
        }
    }
    fn ps(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { .. } => Err(MyErr::new(&format!("no ps in single phase"))),
            Phase::SATVAP { rhog } => Ok(self.alpha.calc(ThermoProp::P, self.T, rhog)),
            Phase::SATLIQ { rhol } => Ok(self.alpha.calc(ThermoProp::P, self.T, rhol)),
            Phase::DOUBLE { x, rhog, rhol } => Ok(x * self.alpha.calc(ThermoProp::P, self.T, rhog)
                + (1.0 - x) * self.alpha.calc(ThermoProp::P, self.T, rhol)),
        }
    }
    fn rhogs(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { .. } => Err(MyErr::new(&format!("no rhogs in single phase"))),
            Phase::SATVAP { rhog } => Ok(rhog),
            Phase::SATLIQ { .. } => Err(MyErr::new(&format!("no rhogs in saturated liquid"))),
            Phase::DOUBLE { rhog, .. } => Ok(rhog),
        }
    }
    fn rhols(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { .. } => Err(MyErr::new(&format!("no rhogs in single phase"))),
            Phase::SATVAP { .. } => Err(MyErr::new(&format!("no rhogs in saturated vapor"))),
            Phase::SATLIQ { rhol } => Ok(rhol),
            Phase::DOUBLE { rhol, .. } => Ok(rhol),
        }
    }
}
