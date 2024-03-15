use super::real_helmholtz_equation::RealHelmholtzEquation;
use super::PropPd;
use super::ThermoProp;
use crate::pubtraits::Flash;
use crate::pubtraits::MyErr;
use crate::pubtraits::Prop;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Read;
use std::path::Path;
///
/// HelmholtzPure 状态方程模型
///
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug)]
pub struct HelmholtzPure {
    #[serde(alias = "helmholtz")]
    eos: RealHelmholtzEquation, // 状态方程模型
    omega: f64, // 偏心因子
    #[serde(skip, default = "default_phase")]
    phase: Phase, // 记录相态
    #[serde(skip, default)]
    T: f64, // 记录温度
}
#[derive(Debug)]
enum Phase {
    SINGLE { rho: f64 },                     // 密度
    SATRHO { rhog: f64, rhol: f64 },         // 饱和气相密度 饱和液相密度
    DOUBLE { rhog: f64, rhol: f64, x: f64 }, // 饱和气相密度 饱和液相密度 干度
}
fn default_phase() -> Phase {
    Phase::SINGLE { rho: 1.0 }
}
pub fn read_json(path: &str) -> Result<HelmholtzPure, MyErr> {
    let mut str_json = String::new();
    let _ = match File::open(&Path::new(path)) {
        Ok(file) => file,
        Err(_) => match File::open(&Path::new(env!("CARGO_MANIFEST_DIR")).join("res").join(path)) {
            Ok(file) => file,
            Err(_) => return Err(MyErr::new(&format!("couldn't find {}", path))),
        },
    }
    .read_to_string(&mut str_json);
    match serde_json::from_str(&str_json) {
        Ok(hp) => Ok(hp),
        Err(_) => Err(MyErr::new(&format!("no alpha(HelmholtzPure) in {}", path))),
    }
}
impl Flash for HelmholtzPure {
    #[allow(non_snake_case)]
    fn t_flash(&mut self, Ts: f64) -> Result<(), MyErr> {
        // Reference [AKASAKA_2008]
        // Newton-Raphson algorithm
        let rhoc = self.eos.Dc();
        let mut rhog: f64;
        let mut rhol: f64;
        let mut Jg: f64;
        let mut Jl: f64;
        let mut Kg: f64;
        let mut Kl: f64;
        let mut DJg: f64;
        let mut DJl: f64;
        let mut DKg: f64;
        let mut DKl: f64;
        let mut Delta: f64;
        let mut is_conv: bool;
        for factor in 0..10 {
            rhog = self.eos.rhogs(Ts) * (1.0 - (factor as f64) / 10.0);
            rhol = self.eos.rhols(Ts) * (1.0 + (factor as f64) / 10.0);
            loop {
                if rhog.is_nan()  // 检查异常结果 NAN
                    || rhol.is_nan()  // 检查异常结果 NAN
                    || self.eos.calc_pd(PropPd::Prho, Ts, rhog) <= 0.0
                    || self.eos.calc_pd(PropPd::Prho, Ts, rhol) <= 0.0
                {
                    is_conv = false;
                    break; // 计算发散 结束循环
                }
                Jg = self.eos.calc_pd(PropPd::J, Ts, rhog);
                Jl = self.eos.calc_pd(PropPd::J, Ts, rhol);
                Kg = self.eos.calc_pd(PropPd::K, Ts, rhog);
                Kl = self.eos.calc_pd(PropPd::K, Ts, rhol);
                if ((Kg - Kl).abs() + (Jg - Jl).abs()) < 1E-8 {
                    // 收敛判据 1E-8 参考[AKASAKA_2008]
                    is_conv = true;
                    break; // 计算收敛 结束循环
                }
                DJg = self.eos.calc_pd(PropPd::Jdelta, Ts, rhog);
                DJl = self.eos.calc_pd(PropPd::Jdelta, Ts, rhol);
                DKg = self.eos.calc_pd(PropPd::Kdelta, Ts, rhog);
                DKl = self.eos.calc_pd(PropPd::Kdelta, Ts, rhol);
                Delta = DJg * DKl - DJl * DKg;
                rhog += ((Kg - Kl) * DJl - (Jg - Jl) * DKl) / Delta * rhoc;
                rhol += ((Kg - Kl) * DJg - (Jg - Jl) * DKg) / Delta * rhoc;
            }
            if is_conv {
                self.phase = Phase::SATRHO { rhog, rhol };
                self.T = Ts;
                return Ok(());
            }
        }
        Err(MyErr::new(&format!("t({})_flash failed!", Ts)))
    }
    #[allow(non_snake_case)]
    fn td_flash(&mut self, T: f64, rho: f64) -> Result<(), MyErr> {
        // 检查是否位于单相区
        if T >= self.eos.Tc() || rho <= 0.85 * self.eos.rhogs(T) || rho >= 1.05 * self.eos.rhols(T)
        {
            self.phase = Phase::SINGLE { rho };
            self.T = T;
            return Ok(());
        }
        // 计算饱和线 判断是否位于两相区
        if let Err(why) = self.t_flash(T) {
            return Err(why);
        } else {
            match self.phase {
                Phase::SATRHO { rhog, rhol } => {
                    if rho < rhog || rho > rhol {
                        // 单相区
                        self.phase = Phase::SINGLE { rho };
                        // self.T = T;  // 温度在饱和线计算过程时就已经设置 所以可以省略
                        return Ok(());
                    } else {
                        // 两相区
                        self.phase = Phase::DOUBLE {
                            rhog,
                            rhol,
                            x: (1.0 / rho - 1.0 / rhol) / (1.0 / rhog + 1.0 / rhol),
                        };
                        // self.T = T;  // 温度在饱和线计算过程时就已经设置 所以可以省略
                        return Ok(());
                    }
                }
                _ => return Err(MyErr::new(&format!("shit-code int td_flash."))),
            }
        }
    }
}
impl Prop for HelmholtzPure {
    fn T(&self) -> Result<f64, MyErr> {
        Ok(self.T)
    }
    fn rho(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(rho),
            Phase::DOUBLE { rhog, rhol, x } => Ok(1.0 / (x / rhog + (1.0 - x) / rhol)),
            Phase::SATRHO { .. } => Err(MyErr::new(&format!("no rho in saturation line"))),
        }
    }
    fn p(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.eos.calc(ThermoProp::P, self.T, rho)),
            Phase::DOUBLE { rhog, rhol, x } => Ok(x * self.eos.calc(ThermoProp::P, self.T, rhog)
                + (1.0 - x) * self.eos.calc(ThermoProp::P, self.T, rhol)),
            Phase::SATRHO { rhog, rhol } => Ok(self.eos.calc(ThermoProp::P, self.T, rhog) / 2.0
                + self.eos.calc(ThermoProp::P, self.T, rhol) / 2.0),
        }
    }
    fn cv(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.eos.calc(ThermoProp::CV, self.T, rho)),
            Phase::DOUBLE { rhog, rhol, x } => Ok(x * self.eos.calc(ThermoProp::CV, self.T, rhog)
                + (1.0 - x) * self.eos.calc(ThermoProp::CV, self.T, rhol)),
            Phase::SATRHO { .. } => Err(MyErr::new(&format!("no cv in saturation line"))),
        }
    }
    fn cp(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.eos.calc(ThermoProp::CP, self.T, rho)),
            Phase::DOUBLE { rhog, rhol, x } => Ok(x * self.eos.calc(ThermoProp::CP, self.T, rhog)
                + (1.0 - x) * self.eos.calc(ThermoProp::P, self.T, rhol)),
            Phase::SATRHO { .. } => Err(MyErr::new(&format!("no cp in saturation line"))),
        }
    }
    fn w(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { rho } => Ok(self.eos.calc(ThermoProp::W, self.T, rho)),
            _ => Err(MyErr::new(&format!("no speed of sound in double phase"))),
        }
    }
    fn ps(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { .. } => Err(MyErr::new(&format!("no ps in single phase"))),
            Phase::DOUBLE { rhog, rhol, x } => Ok(x * self.eos.calc(ThermoProp::P, self.T, rhog)
                + (1.0 - x) * self.eos.calc(ThermoProp::P, self.T, rhol)),
            Phase::SATRHO { rhog, rhol } => Ok(self.eos.calc(ThermoProp::P, self.T, rhog) / 2.0
                + self.eos.calc(ThermoProp::P, self.T, rhol) / 2.0),
        }
    }
    fn rhogs(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { .. } => Err(MyErr::new(&format!("no rhogs in single phase"))),
            Phase::SATRHO { rhog, .. } => Ok(rhog),
            Phase::DOUBLE { rhog, .. } => Ok(rhog),
        }
    }
    fn rhols(&self) -> Result<f64, MyErr> {
        match self.phase {
            Phase::SINGLE { .. } => Err(MyErr::new(&format!("no rhogs in single phase"))),
            Phase::SATRHO { rhol, .. } => Ok(rhol),
            Phase::DOUBLE { rhol, .. } => Ok(rhol),
        }
    }
}
// 单元测试
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test() {
        let mut so2: HelmholtzPure = read_json("SO2.json").expect("no SO2");
        let serialized = serde_json::to_string_pretty(&so2).unwrap();
        println!("serialized = {}", serialized);
        println!("read_json pass!");
        for t in 197..432 {
            match so2.t_flash(t as f64) {
                Ok(()) => println!("Ts={}\tPs={}", t, so2.ps().unwrap()),
                Err(_) => {
                    println!("error occured in {}K.", t);
                    break;
                }
            }
        }
        println!("t_flash pass!");
        match so2.td_flash(350.0, 10000.0) {
            Ok(()) => println!("PS={}", so2.ps().expect("no ps")),
            Err(_) => {
                println!("error occured in t_flash.");
            }
        }
        println!("td_flash pass!");
    }
}
