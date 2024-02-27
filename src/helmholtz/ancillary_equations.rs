use serde::{Deserialize, Serialize};
///
/// 辅助方程：
/// 饱和蒸汽压方程
/// 饱和气相密度方程
/// 饱和液相密度方程
///
#[derive(Serialize, Deserialize, Debug)]
struct AuxiliaryEquationTerm {
    n: f64, // 系数
    t: f64, // 指数
}
impl AuxiliaryEquationTerm {
    fn calc(&self, theta: f64) -> f64 {
        self.n * theta.powf(self.t)
    }
}
#[derive(Serialize, Deserialize, Debug)]
pub struct VaporPressureEquation {
    flag: i32, // 方程标记
    terms: Vec<AuxiliaryEquationTerm>,
}
impl VaporPressureEquation {
    #[allow(non_snake_case)]
    pub fn calc(&self, T: f64, Tr: f64, pr: f64) -> f64 {
        let theta = T / Tr;
        let mut p_pr = 0.0;
        pr * (match self.flag {
            1 => {
                for i in self.terms.iter() {
                    p_pr += i.calc(theta);
                }
                (p_pr / theta).exp()
            }
            _ => {
                println!("no flag={} in vapor_pressure_equation\n", self.flag);
                1.0
            }
        })
    }
}
#[derive(Serialize, Deserialize, Debug)]
pub struct SaturatedVaporDensityEquation {
    flag: i32, // 方程标记
    terms: Vec<AuxiliaryEquationTerm>,
}
impl SaturatedVaporDensityEquation {
    #[allow(non_snake_case)]
    pub fn calc(&self, T: f64, Tr: f64, rhor: f64) -> f64 {
        let theta = T / Tr;
        let mut rho_rhor = 0.0;
        rhor * (match self.flag {
            1 => {
                for i in self.terms.iter() {
                    rho_rhor += i.calc(theta);
                }
                rho_rhor.exp()
            }
            2 => {
                for i in self.terms.iter() {
                    rho_rhor += i.calc(theta / 3.0);
                }
                rho_rhor.exp()
            }
            _ => {
                println!(
                    "no flag={} in saturated_vapor_density_equation\n",
                    self.flag
                );
                1.0
            }
        })
    }
}
#[derive(Serialize, Deserialize, Debug)]
pub struct SaturatedLiquidDensityEquation {
    flag: i32, // 方程标记
    terms: Vec<AuxiliaryEquationTerm>,
}
impl SaturatedLiquidDensityEquation {
    #[allow(non_snake_case)]
    pub fn calc(&self, T: f64, Tr: f64, rhor: f64) -> f64 {
        let theta = T / Tr;
        let mut rho_rhor = 0.0;
        rhor * (match self.flag {
            1 => {
                for i in self.terms.iter() {
                    rho_rhor += i.calc(theta);
                }
                rho_rhor + 1.0
            }
            2 => {
                for i in self.terms.iter() {
                    rho_rhor += i.calc(theta / 3.0);
                }
                rho_rhor + 1.0
            }
            3 => {
                for i in self.terms.iter() {
                    rho_rhor += i.calc(theta / 3.0);
                }
                rho_rhor.exp()
            }
            _ => {
                println!(
                    "no flag={} in saturated_liquid_density_equation\n",
                    self.flag
                );
                1.0
            }
        })
    }
}
