use anyhow::anyhow;
use lazy_static::lazy_static;
use pyo3::{pyclass, pymethods};
use std::collections::HashMap;
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
use eta::METALS_TO_ETAPARAMS;
use lambda::METALS_TO_LAMBDAPARAMS;
use rho::METALS_TO_RHOPARAMS;
/// liquid metal
/// + Density, unit: kg/m3
/// + Thermal Conductivity, unit: W/m/K
/// + Viscosity, unit: mPa*s
#[pyclass]
#[allow(non_snake_case)]
pub struct LiquidMetal {
    metal: Metals,
}
#[pymethods]
#[allow(non_snake_case)]
impl LiquidMetal {
    #[new]
    pub fn new_metal(name: &str) -> anyhow::Result<LiquidMetal> {
        if let Some(metal) = STRING_TO_METALS.get(name) {
            Ok(Self {
                metal: metal.clone(),
            })
        } else {
            Err(anyhow!(LiquidMetalErr::NoLiquidMetal))
        }
    }
    /// calculate density of liquid metals at 0.1 MPa, unit: kg/m3
    pub fn calc_rho(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(rho_params) = METALS_TO_RHOPARAMS.get(&self.metal) {
            rho_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
    /// calculate thermal conductivity of liquid metals at 0.1 MPa, unit: W/m/K
    pub fn calc_lambda(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(lambda_params) = METALS_TO_LAMBDAPARAMS.get(&self.metal) {
            lambda_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
    /// calculate viscosity of liquid metals at 0.1 MPa, unit: mPa*s
    pub fn calc_eta(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(eta_params) = METALS_TO_ETAPARAMS.get(&self.metal) {
            eta_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
}
lazy_static! {
    static ref STRING_TO_METALS: HashMap<String, Metals> = HashMap::from([
        (String::from("Al"), Metals::Al),
        (String::from("Si"), Metals::Si),
        (String::from("Ti"), Metals::Ti),
        (String::from("V"), Metals::V),
        (String::from("Cr"), Metals::Cr),
        (String::from("Fe"), Metals::Fe),
        (String::from("Co"), Metals::Co),
        (String::from("Ni"), Metals::Ni),
        (String::from("Cu"), Metals::Cu),
        (String::from("Zn"), Metals::Zn),
        (String::from("Ga"), Metals::Ga),
        (String::from("Ge"), Metals::Ge),
        (String::from("Zr"), Metals::Zr),
        (String::from("Nb"), Metals::Nb),
        (String::from("Mo"), Metals::Mo),
        (String::from("Ag"), Metals::Ag),
        (String::from("Cd"), Metals::Cd),
        (String::from("In"), Metals::In),
        (String::from("Sn"), Metals::Sn),
        (String::from("Sb"), Metals::Sb),
        (String::from("Hf"), Metals::Hf),
        (String::from("Ta"), Metals::Ta),
        (String::from("W"), Metals::W),
        (String::from("Hg"), Metals::Hg),
        (String::from("Tl"), Metals::Tl),
        (String::from("Pb"), Metals::Pb),
        (String::from("Bi"), Metals::Bi),
    ]);
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
