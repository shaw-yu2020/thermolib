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
/// + Density
/// + Thermal Conductivity
/// + Viscosity
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
    pub fn full_name(&self) -> anyhow::Result<String> {
        if let Some(metal) = METALS_TO_STRING.get(&self.metal) {
            Ok(metal.clone())
        } else {
            Err(anyhow!(LiquidMetalErr::NoLiquidMetal))
        }
    }
    /// calculate density of liquid metals at 0.1 MPa
    pub fn calc_rho(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(rho_params) = METALS_TO_RHOPARAMS.get(&self.metal) {
            rho_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
    /// calculate thermal conductivity of liquid metals at 0.1 MPa
    pub fn calc_lambda(&self, T: f64) -> anyhow::Result<f64> {
        if let Some(lambda_params) = METALS_TO_LAMBDAPARAMS.get(&self.metal) {
            lambda_params.calc(T)
        } else {
            Err(anyhow!(LiquidMetalErr::NoProperty))
        }
    }
    /// calculate viscosity of liquid metals at 0.1 MPa
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
        (String::from("Ti"), Metals::Ti),
        (String::from("V"), Metals::V),
        (String::from("Cr"), Metals::Cr),
        (String::from("Fe"), Metals::Fe),
        (String::from("Zr"), Metals::Zr),
        (String::from("Nb"), Metals::Nb),
        (String::from("Mo"), Metals::Mo),
        (String::from("Hf"), Metals::Hf),
        (String::from("Ta"), Metals::Ta),
        (String::from("W"), Metals::W),
    ]);
}
#[derive(Clone, Hash, Eq, PartialEq)]
pub enum Metals {
    Al, // 13
    Ti, // 22
    V,  // 23
    Cr, // 24
    Fe, // 26
    Zr, // 40
    Nb, // 41
    Mo, // 42
    Hf, // 72
    Ta, // 73
    W,  // 74
}
lazy_static! {
    static ref METALS_TO_STRING: HashMap<Metals, String> = HashMap::from([
        (Metals::Al, String::from("Aluminum")),
        (Metals::Ti, String::from("Titanium")),
        (Metals::V, String::from("Vanadium")),
        (Metals::Cr, String::from("Chromium")),
        (Metals::Fe, String::from("Iron")),
        (Metals::Zr, String::from("Zirconium")),
        (Metals::Nb, String::from("Niobium")),
        (Metals::Mo, String::from("Molybdenum")),
        (Metals::Hf, String::from("Hafnium")),
        (Metals::Ta, String::from("Tantalum")),
        (Metals::W, String::from("Tungsten"))
    ]);
}
