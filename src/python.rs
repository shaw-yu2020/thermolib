/// A Python module implemented in Rust.
use pyo3::prelude::*;
#[pyfunction]
fn hello() -> PyResult<String> {
    Ok(String::from("Hello, Rust And Python."))
}
use crate::Helmholtz;
use crate::IdealGas;
use crate::LiquidMetal;
use crate::PcSaftGlyPure;
use crate::PcSaftPure;
use crate::Pr;
use crate::Rk;
use crate::Srk;
use crate::Vdw;
/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "thermolib")] // Rename pymodule
pub fn pylib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hello, m)?)?;
    m.add_class::<Vdw>()?;
    m.add_class::<Rk>()?;
    m.add_class::<Srk>()?;
    m.add_class::<Pr>()?;
    m.add_class::<Helmholtz>()?;
    m.add_class::<LiquidMetal>()?;
    m.add_class::<PcSaftGlyPure>()?;
    m.add_class::<PcSaftPure>()?;
    m.add_class::<IdealGas>()?;
    Ok(())
}
