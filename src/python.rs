use crate::pubtraits::Flash;
use crate::pubtraits::Prop;
use crate::HelmholtzPure;
use crate::Pr;
use pyo3::prelude::*;
/// A Python module implemented in Rust.
/// The name of this function must match the `lib.name` in the `Cargo.toml`,
/// else Python will not be able to import the module.
#[pymodule] // PYTHON 模块
#[pyo3(name = "thermolib")] // 调整 名称
pub fn pylib(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<NewModel>()?;
    Ok(())
}
#[pyclass] // PYTHON 类
pub struct NewModel {
    model: Models,
}
enum Models {
    PR(Pr),
    HELMHOLTZ(HelmholtzPure),
}
#[pymethods] // PYTHON 方法
impl NewModel {
    #[new]
    fn new(model: &str) -> Self {
        match model {
            "PR" => NewModel {
                model: (Models::PR(Pr::new_fluid(430.64, 7.884E6, 8.314462618, 0.2557))),
            },
            "HELMHOLTZ" => NewModel {
                model: (Models::HELMHOLTZ(
                    HelmholtzPure::read_json("SO2.json").expect("no SO2.json"),
                )),
            },
            _ => {
                println!("error input, use default = PR(SO2).");
                NewModel {
                    model: (Models::PR(Pr::new_fluid(430.64, 7.884E6, 8.314462618, 0.2557))),
                }
            }
        }
    }
    #[allow(non_snake_case)]
    fn t_flash(&mut self, T: f64) -> PyResult<()> {
        Ok(match &mut self.model {
            Models::PR(model) => model.t_flash(T).expect("REASON"),
            Models::HELMHOLTZ(model) => model.t_flash(T).expect("REASON"),
        })
    }
    #[allow(non_snake_case)]
    fn td_flash(&mut self, T: f64, rho: f64) -> PyResult<()> {
        Ok(match &mut self.model {
            Models::PR(model) => model.td_flash(T, rho).expect("REASON"),
            Models::HELMHOLTZ(model) => model.td_flash(T, rho).expect("REASON"),
        })
    }
    #[allow(non_snake_case)]
    fn tp_flash(&mut self, T: f64, p: f64) -> PyResult<()> {
        Ok(match &mut self.model {
            Models::PR(model) => model.tp_flash(T, p).expect("REASON"),
            Models::HELMHOLTZ(model) => model.tp_flash(T, p).expect("REASON"),
        })
    }
    #[allow(non_snake_case)]
    fn T(&self) -> f64 {
        match &self.model {
            Models::PR(model) => model.T().expect("REASON"),
            Models::HELMHOLTZ(model) => model.T().expect("REASON"),
        }
    }
    fn rho(&self) -> f64 {
        match &self.model {
            Models::PR(model) => model.rho().expect("REASON"),
            Models::HELMHOLTZ(model) => model.rho().expect("REASON"),
        }
    }
    fn p(&self) -> f64 {
        match &self.model {
            Models::PR(model) => model.p().expect("REASON"),
            Models::HELMHOLTZ(model) => model.p().expect("REASON"),
        }
    }
}
