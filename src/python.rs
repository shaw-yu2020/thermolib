/// A Python module implemented in Rust.
use pyo3::prelude::*;
#[pyfunction]
fn hello() -> PyResult<String> {
    Ok(String::from("Hello, Rust And Python."))
}
use crate::PcSaftGlyPure;
/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "thermolib")] // Rename pymodule
fn pylib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hello, m)?)?;
    m.add_class::<PcSaftGlyPure>()?;
    Ok(())
}
