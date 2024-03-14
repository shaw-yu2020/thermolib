use super::ancillary_equations::SaturatedLiquidDensityEquation;
use super::ancillary_equations::SaturatedVaporDensityEquation;
use super::ancillary_equations::VaporPressureEquation;
use super::real_helmholtz_equation::RealHelmholtzEquation;
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
    alpha: RealHelmholtzEquation,
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
        Ok(_) => {
            print!("{} contains:\n{}\n", path.display(), str_json);
            println!("{} is loaded.", path.display())
        }
        Err(why) => panic!("couldn't read {}: {:?}", path.display(), why),
    }
    serde_json::from_str(&str_json).unwrap()
}
