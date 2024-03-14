use super::ancillary_equations::SaturatedLiquidDensityEquation;
use super::ancillary_equations::SaturatedVaporDensityEquation;
use super::ancillary_equations::VaporPressureEquation;
use super::real_helmholtz_equation::RealHelmholtzEquation;
use super::ThermoProp;
use crate::pubtraits::MyErr;
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
