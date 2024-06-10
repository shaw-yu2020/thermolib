/// A Python module implemented in Rust.
use pyo3::prelude::*;
#[pyfunction]
fn hello() -> PyResult<String> {
    Ok(String::from("Hello, Rust And Python."))
}
use crate::Helmholtz;
use crate::Pr;
use crate::Rk;
use crate::Srk;
use crate::Vdw;
/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "thermolib")] // Rename pymodule
pub fn pylib(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hello, m)?)?;
    m.add_class::<Vdw>()?;
    m.add_class::<Rk>()?;
    m.add_class::<Srk>()?;
    m.add_class::<Pr>()?;
    m.add_class::<Helmholtz>()?;
    Ok(())
}
