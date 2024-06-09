use serde::Deserialize;
/// Ancillary Equations
/// + Vapor Pressure Equation
/// + Saturated Vapor Density Equation
/// + Saturated Liquid Density Equation
#[derive(Deserialize)]
struct AncEqnTerm {
    n: f64, // coefficient
    t: f64, // exponent
}
impl AncEqnTerm {
    fn calc(&self, theta: f64) -> f64 {
        self.n * theta.powf(self.t)
    }
}
/// Vapor Pressure Equation
#[derive(Deserialize)]
pub struct PsEqn {
    flag: i32,
    terms: Vec<AncEqnTerm>,
}
impl PsEqn {
    #[allow(non_snake_case)]
    pub fn calc(&self, T: f64, Tr: f64, pr: f64) -> f64 {
        let theta = 1.0 - T / Tr;
        pr * (match self.flag {
            1 => (self.terms.iter().map(|term| term.calc(theta)).sum::<f64>() * Tr / T).exp(),
            _ => {
                println!("no flag={} in vapor_pressure_equation", self.flag);
                1.0
            }
        })
    }
}
/// Saturated Vapor Density Equation
#[derive(Deserialize)]
pub struct RhovEqn {
    flag: i32,
    terms: Vec<AncEqnTerm>,
}
impl RhovEqn {
    #[allow(non_snake_case)]
    pub fn calc(&self, T: f64, Tr: f64, rhor: f64) -> f64 {
        let theta = 1.0 - T / Tr;
        rhor * (match self.flag {
            1 => self
                .terms
                .iter()
                .map(|term| term.calc(theta))
                .sum::<f64>()
                .exp(),
            _ => {
                println!("no flag={} in saturated_vapor_density_equation", self.flag);
                1.0
            }
        })
    }
}
/// Saturated Liquid Density Equation
#[derive(Deserialize)]
pub struct RholEqn {
    flag: i32,
    terms: Vec<AncEqnTerm>,
}
impl RholEqn {
    #[allow(non_snake_case)]
    pub fn calc(&self, T: f64, Tr: f64, rhor: f64) -> f64 {
        let theta = 1.0 - T / Tr;
        rhor * (match self.flag {
            1 => self.terms.iter().map(|term| term.calc(theta)).sum::<f64>() + 1.0,
            _ => {
                println!("no flag={} in saturated_liquid_density_equation", self.flag);
                1.0
            }
        })
    }
}
