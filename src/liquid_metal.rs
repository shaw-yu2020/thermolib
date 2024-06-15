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
        if let Some(metal) = MAP_SM.get(name) {
            Ok(Self {
                metal: metal.clone(),
            })
        } else {
            Err(anyhow!(LiquidMetalErr::NoLiquidMetal))
        }
    }
    pub fn full_name(&self) -> anyhow::Result<String> {
        if let Some(metal) = MAP_MS.get(&self.metal) {
            Ok(metal.clone())
        } else {
            Err(anyhow!(LiquidMetalErr::NoLiquidMetal))
        }
    }
    pub fn calc_rho(&self, T: f64) -> anyhow::Result<f64> {
        self.rho(T)
    }
    pub fn calc_lambda(&self, T: f64) -> anyhow::Result<f64> {
        self.lambda(T)
    }
    pub fn calc_eta(&self, T: f64) -> anyhow::Result<f64> {
        self.eta(T)
    }
}
lazy_static! {
    static ref MAP_SM: HashMap<String, Metals> = HashMap::from([
        (String::from("Ti"), Metals::Ti),
        (String::from("V"), Metals::V),
        (String::from("Cr"), Metals::Cr),
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
    Ti, // 22
    V,  // 23
    Cr, // 24
    Zr, // 40
    Nb, // 41
    Mo, // 42
    Hf, // 72
    Ta, // 73
    W,  // 74
}
lazy_static! {
    static ref MAP_MS: HashMap<Metals, String> = HashMap::from([
        (Metals::Ti, String::from("Titanium")),
        (Metals::V, String::from("Vanadium")),
        (Metals::Cr, String::from("Chromium")),
        (Metals::Zr, String::from("Zirconium")),
        (Metals::Nb, String::from("Niobium")),
        (Metals::Mo, String::from("Molybdenum")),
        (Metals::Hf, String::from("Hafnium")),
        (Metals::Ta, String::from("Tantalum")),
        (Metals::W, String::from("Tungsten"))
    ]);
}
