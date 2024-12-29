use anyhow::anyhow;
#[cfg(feature = "with_pyo3")]
use pyo3::{pyclass, pymethods};
use std::collections::HashMap;
use std::sync::OnceLock;
use thiserror::Error;
#[derive(Debug, Error)]
enum LiquidMetalErr {
    #[error("no liquid metal")]
    NoLiquidMetal,
    #[error("no property")]
    NoProperty,
    #[error("T is too min")]
    TisTooMin,
    #[error("T is too max")]
    TisTooMax,
}
mod eta;
mod lambda;
mod rho;
use eta::metals2etaparams;
use lambda::metals2lambdaparams;
use rho::metals2rhoparams;
/// liquid metal
/// + Density, unit: kg/m3
/// + Thermal Conductivity, unit: W/m/K
/// + Viscosity, unit: mPa*s
/// ```
/// use thermolib::LiquidMetal;
/// let T = 1800.0; // K
/// let metal = LiquidMetal::new_metal("Si").unwrap();
/// let rho = metal.calc_rho(T).unwrap(); // kg/m3
/// assert_eq!(2520.0, (rho * 1e0).round() / 1e0);
/// let eta = metal.calc_eta(T).unwrap(); // mPa*s
/// assert_eq!(0.541, (eta * 1e3).round() / 1e3);
/// let lambda = metal.calc_lambda(T).unwrap(); // W/m/K
/// assert_eq!(54.88, (lambda * 1e2).round() / 1e2);
/// ```
#[cfg_attr(feature = "with_pyo3", pyclass)]
#[allow(non_snake_case)]
pub struct LiquidMetal {
    metal: Metals,
}
#[allow(non_snake_case)]
impl LiquidMetal {
    pub fn new_metal(name: &str) -> anyhow::Result<LiquidMetal> {
        if let Some(metal) = string2metals().get(name) {
            Ok(Self {
                metal: metal.clone(),
            })
        } else {
            Err(anyhow!(LiquidMetalErr::NoLiquidMetal))
        }
    }
}
#[cfg_attr(feature = "with_pyo3", pymethods)]
#[allow(non_snake_case)]
impl LiquidMetal {
    #[cfg(feature = "with_pyo3")]
    #[new]
    pub fn new_py(name: &str) -> anyhow::Result<LiquidMetal> {
        Self::new_metal(name)
    }
    /// calculate density of liquid metals at 0.1 MPa, unit: kg/m3
    pub fn calc_rho(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(rho_params) = metals2rhoparams().get(&self.metal) {
            rho_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
    /// calculate thermal conductivity of liquid metals at 0.1 MPa, unit: W/m/K
    pub fn calc_lambda(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(lambda_params) = metals2lambdaparams().get(&self.metal) {
            lambda_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
    /// calculate viscosity of liquid metals at 0.1 MPa, unit: mPa*s
    pub fn calc_eta(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(eta_params) = metals2etaparams().get(&self.metal) {
            eta_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
}
fn string2metals() -> &'static HashMap<String, Metals> {
    static STRING_TO_METALS: OnceLock<HashMap<String, Metals>> = OnceLock::new();
    STRING_TO_METALS.get_or_init(|| {
        let mut hm = HashMap::new();
        hm.insert(String::from("Al"), Metals::Al);
        hm.insert(String::from("Si"), Metals::Si);
        hm.insert(String::from("Ti"), Metals::Ti);
        hm.insert(String::from("V"), Metals::V);
        hm.insert(String::from("Cr"), Metals::Cr);
        hm.insert(String::from("Fe"), Metals::Fe);
        hm.insert(String::from("Co"), Metals::Co);
        hm.insert(String::from("Ni"), Metals::Ni);
        hm.insert(String::from("Cu"), Metals::Cu);
        hm.insert(String::from("Zn"), Metals::Zn);
        hm.insert(String::from("Ga"), Metals::Ga);
        hm.insert(String::from("Ge"), Metals::Ge);
        hm.insert(String::from("Zr"), Metals::Zr);
        hm.insert(String::from("Nb"), Metals::Nb);
        hm.insert(String::from("Mo"), Metals::Mo);
        hm.insert(String::from("Ag"), Metals::Ag);
        hm.insert(String::from("Cd"), Metals::Cd);
        hm.insert(String::from("In"), Metals::In);
        hm.insert(String::from("Sn"), Metals::Sn);
        hm.insert(String::from("Sb"), Metals::Sb);
        hm.insert(String::from("Hf"), Metals::Hf);
        hm.insert(String::from("Ta"), Metals::Ta);
        hm.insert(String::from("W"), Metals::W);
        hm.insert(String::from("Hg"), Metals::Hg);
        hm.insert(String::from("Tl"), Metals::Tl);
        hm.insert(String::from("Pb"), Metals::Pb);
        hm.insert(String::from("Bi"), Metals::Bi);
        hm
    })
}
#[derive(Clone, Hash, Eq, PartialEq)]
pub enum Metals {
    Al, // 13
    Si, // 14
    Ti, // 22
    V,  // 23
    Cr, // 24
    Fe, // 26
    Co, // 27
    Ni, // 28
    Cu, // 29
    Zn, // 30
    Ga, // 31
    Ge, // 32
    Zr, // 40
    Nb, // 41
    Mo, // 42
    Ag, // 47
    Cd, // 48
    In, // 49
    Sn, // 50
    Sb, // 51
    Hf, // 72
    Ta, // 73
    W,  // 74
    Hg, // 80
    Tl, // 81
    Pb, // 82
    Bi, // 83
}
