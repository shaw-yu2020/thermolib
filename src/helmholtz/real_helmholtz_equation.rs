use super::ancillary_equations::SaturatedLiquidDensityEquation;
use super::ancillary_equations::SaturatedVaporDensityEquation;
use super::ancillary_equations::VaporPressureEquation;
use super::ideal_helmholtz_equation::IdealHelmholtzEquation;
use super::residual_helmholtz_equation::ResidualHelmholtzEquation;
use super::Alpha0Dtau;
use super::AlpharDD;
use super::PropPd;
use super::ThermoProp;
use serde::{Deserialize, Serialize};
///
/// 实际气体亥姆霍兹方程
///
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug)]
pub struct RealHelmholtzEquation {
    Tc: f64, // 临界温度
    Pc: f64, // 临界压力
    Dc: f64, // 临界密度
    R: f64,  // 气体常数
    M: f64,  // 摩尔质量
    alpha0: IdealHelmholtzEquation,
    alphar: ResidualHelmholtzEquation,
    ps: VaporPressureEquation,
    rhogs: SaturatedVaporDensityEquation,
    rhols: SaturatedLiquidDensityEquation,
}
impl RealHelmholtzEquation {
    #[allow(non_snake_case)]
    pub fn calc(&self, tp: ThermoProp, T: f64, D: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = D / self.Dc;
        match tp {
            ThermoProp::Z => 1.0 + self.alphar.calc(AlpharDD::D01, tau, delta),
            ThermoProp::P => (D * self.R * T) * (1.0 + self.alphar.calc(AlpharDD::D01, tau, delta)),
            ThermoProp::CV => {
                self.R
                    * (-self.alpha0.calc(Alpha0Dtau::D2, tau, self.Tc)
                        - self.alphar.calc(AlpharDD::D20, tau, delta))
            }
            ThermoProp::CP => {
                self.R
                    * (-self.alpha0.calc(Alpha0Dtau::D2, tau, self.Tc)
                        - self.alphar.calc(AlpharDD::D20, tau, delta)
                        + (1.0 + self.alphar.calc(AlpharDD::D01, tau, delta)
                            - self.alphar.calc(AlpharDD::D11, tau, delta))
                        .powi(2)
                            / (1.0
                                + 2.0 * self.alphar.calc(AlpharDD::D01, tau, delta)
                                + self.alphar.calc(AlpharDD::D02, tau, delta)))
            }
            ThermoProp::W => {
                let Rg = if self.R < 10.0 {
                    self.R / self.M
                } else {
                    self.R
                };
                (Rg * T
                    * ((1.0
                        + 2.0 * self.alphar.calc(AlpharDD::D01, tau, delta)
                        + self.alphar.calc(AlpharDD::D02, tau, delta))
                        - (1.0 + self.alphar.calc(AlpharDD::D01, tau, delta)
                            - self.alphar.calc(AlpharDD::D11, tau, delta))
                        .powi(2)
                            / (self.alpha0.calc(Alpha0Dtau::D2, tau, self.Tc)
                                + self.alphar.calc(AlpharDD::D20, tau, delta))))
                .sqrt()
            }
            ThermoProp::S => {
                self.R
                    * (tau * self.alpha0.calc(Alpha0Dtau::D1, tau, self.Tc)
                        + self.alphar.calc(AlpharDD::D10, tau, delta)
                        - (self.alpha0.calc(Alpha0Dtau::D0, tau, self.Tc) + delta.ln())
                        - self.alphar.calc(AlpharDD::D00, tau, delta))
            }
            ThermoProp::U => {
                (self.R * T)
                    * (self.alpha0.calc(Alpha0Dtau::D1, tau, self.Tc)
                        + self.alphar.calc(AlpharDD::D10, tau, delta))
            }
            ThermoProp::H => {
                (self.R * T)
                    * (1.0
                        + self.alpha0.calc(Alpha0Dtau::D1, tau, self.Tc)
                        + self.alphar.calc(AlpharDD::D10, tau, delta)
                        + self.alphar.calc(AlpharDD::D01, tau, delta))
            }
            ThermoProp::A => {
                (self.R * T)
                    * ((self.alpha0.calc(Alpha0Dtau::D0, tau, self.Tc) + delta.ln())
                        + self.alphar.calc(AlpharDD::D00, tau, delta))
            }
            ThermoProp::G => {
                (self.R * T)
                    * (1.0
                        + (self.alpha0.calc(Alpha0Dtau::D0, tau, self.Tc) + delta.ln())
                        + self.alphar.calc(AlpharDD::D00, tau, delta)
                        + self.alphar.calc(AlpharDD::D01, tau, delta))
            }
        }
    }
    #[allow(non_snake_case)]
    pub fn calc_pd(&self, pd: PropPd, T: f64, D: f64) -> f64 {
        let tau = self.Tc / T;
        let delta = D / self.Dc;
        match pd {
            PropPd::J => delta * (1.0 + self.alphar.calc(AlpharDD::D01, tau, delta)),
            PropPd::K => {
                self.alphar.calc(AlpharDD::D01, tau, delta)
                    + self.alphar.calc(AlpharDD::D00, tau, delta)
                    + delta.ln()
            }
            PropPd::Jdelta => {
                1.0 + 2.0 * self.alphar.calc(AlpharDD::D01, tau, delta)
                    + self.alphar.calc(AlpharDD::D02, tau, delta)
            }
            PropPd::Kdelta => {
                2.0 * self.alphar.calc(AlpharDD::D01, tau, delta)
                    + self.alphar.calc(AlpharDD::D02, tau, delta)
                    + 1.0 / delta
            }
        }
    }
    #[allow(non_snake_case)]
    pub fn ps(&self, T: f64) -> f64 {
        self.ps.calc(T, self.Tc, self.Pc)
    }
    #[allow(non_snake_case)]
    pub fn rhogs(&self, T: f64) -> f64 {
        self.rhogs.calc(T, self.Tc, self.Dc)
    }
    #[allow(non_snake_case)]
    pub fn rhols(&self, T: f64) -> f64 {
        self.rhols.calc(T, self.Tc, self.Dc)
    }
}
