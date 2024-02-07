use crate::pubtraits::Flash;
use crate::pubtraits::MyErr;
use crate::pubtraits::Prop;
///
/// PR 方程
///
#[allow(non_snake_case)]
pub struct Pr {
    Tc: f64,               // 临界温度 Tc: K
    pc: f64,               // 临界压力 pc: Pa
    R: f64,                // 理想气体常数 R: J/mol/K || J/kg/K
    kappa: f64,            // 根据偏心因子 omega 计算 kappa
    ac: f64,               // 临界点的 a
    bc: f64,               // 临界点的 b
    Zc: f64,               // 临界点的 Z
    T: f64,                // 单相区的温度 T&P 计算 A&B
    P: f64,                // 单相区的压力 P&T 计算 A&B
    A: f64,                // 用于计算 Z
    B: f64,                // 用于计算 Z
    Z: f64,                // 单相区的压缩因子 Z
    Zg: f64,               // 饱和气相压缩因子 Zg
    Zl: f64,               // 饱和液相压缩因子 Zl
    is_single_phase: bool, // 是否为单相区
}
fn zero_fluid() -> Pr {
    Pr {
        Tc: 1.0,    // 外部变量
        pc: 1.0,    // 外部变量
        R: 1.0,     // 外部变量
        kappa: 1.0, // 外部变量
        ac: 1.0,    // 外部变量
        bc: 1.0,    // 外部变量
        Zc: 0.307,  // 内部变量 0.307
        T: 0.0,     // 内部变量 0.0
        P: 0.0,     // 内部变量 0.0
        A: 0.0,     // 内部变量 0.0
        B: 0.0,     // 内部变量 0.0
        Z: 0.0,     // 内部变量 0.0
        Zg: 0.0,    // 内部变量 0.0
        Zl: 0.0,    // 内部变量 0.0
        is_single_phase: true,
    }
}
#[allow(non_snake_case)]
impl Pr {
    pub fn new_fluid(Tc: f64, pc: f64, R: f64, omega: f64) -> Pr {
        let mut fluid = Pr {
            Tc,
            pc,
            R,
            kappa: 0.37464 + 1.54226 * omega - 0.26992 * omega.powi(2),
            ac: 0.45724 * R.powi(2) * Tc.powi(2) / pc,
            bc: 0.07780 * R * Tc / pc,
            ..zero_fluid()
        };
        fluid.tp_flash(273.15, 0.1E6).unwrap();
        fluid
    }
    pub fn set_fluid(&mut self, Tc: f64, pc: f64, R: f64, omega: f64) {
        *self = Pr {
            Tc,
            pc,
            R,
            kappa: 0.37464 + 1.54226 * omega - 0.26992 * omega.powi(2),
            ac: 0.45724 * R.powi(2) * Tc.powi(2) / pc,
            bc: 0.07780 * R * Tc / pc,
            ..zero_fluid()
        };
        self.tp_flash(273.15, 0.1E6).unwrap();
    }
    fn calc_AB(&mut self, T: f64, p: f64) {
        self.T = T;
        self.P = p;
        let a = self.ac * (1.0 + self.kappa * (1.0 - (T / self.Tc).sqrt())).powi(2);
        self.A = a * p / self.R.powi(2) / self.T.powi(2);
        self.B = self.bc * p / self.R / self.T; // b=bc
    }
    fn calc_root(&mut self) -> bool {
        // 根据 A&B 计算 Z||Zg&Zl 返回结果->是否在单相区
        let b = -(1.0 - self.B);
        let c = self.A - 3.0 * self.B.powi(2) - 2.0 * self.B;
        let d = -(self.A * self.B - self.B.powi(2) - self.B.powi(3));
        let A = b.powi(2) - 3.0 * c;
        let B = b * c - 9.0 * d;
        let C = c.powi(2) - 3.0 * b * d;
        let Delta = B.powi(2) - 4.0 * A * C;
        if Delta > 0.0 {
            // Delta>0 方程有一个正实根和一对共轭复根 仅仅保留实根！!
            let Y1 = A * b + 1.5 * (-B + Delta.sqrt());
            let Y2 = A * b + 1.5 * (-B - Delta.sqrt());
            self.Z = (-b - (Y1.cbrt() + Y2.cbrt())) / 3.0;
            return true;
        } else if Delta < 0.0 {
            // Delta<0 方程有三个不相等的实根 舍弃中间根！！
            let theta3 = ((2.0 * A * b - 3.0 * B) / (2.0 * A * A.sqrt())).acos() / 3.0;
            let sqrt_A = f64::sqrt(A);
            let sqrt_3 = f64::sqrt(3.0);
            let x1 = (-b - 2.0 * sqrt_A * theta3.cos()) / 3.0;
            let x2 = (-b + sqrt_A * (theta3.cos() + sqrt_3 * theta3.sin())) / 3.0;
            let x3 = (-b + sqrt_A * (theta3.cos() - sqrt_3 * theta3.sin())) / 3.0;
            let Zg = x1.max(x2).max(x3);
            let Zl = x1.min(x2).min(x3);
            // 避免压缩因子出现负值
            if Zl < 0.0 {
                self.Z = Zg;
                return true;
            } else {
                self.Zg = Zg;
                self.Zl = Zl;
                return false;
            }
        } else {
            // Delta=0 方程有三个正实根 其中有一个二重根 舍弃二重根！！
            self.Z = B / A - b;
            return true;
        }
    }
    fn calc_lnfp(&self, Z: f64) -> f64 {
        let sqrt2 = f64::sqrt(2.0);
        (Z - 1.0)
            - (Z - self.B).ln()
            - self.A / (2.0 * sqrt2 * self.B)
                * ((Z + self.B + sqrt2 * self.B) / (Z + self.B - sqrt2 * self.B)).ln()
    }
    fn calc_diff_lnfpgl(&mut self) -> f64 {
        if !self.calc_root() {
            self.calc_lnfp(self.Zg) - self.calc_lnfp(self.Zl)
        } else {
            if self.Z < self.Zc {
                -self.calc_lnfp(self.Z) // 液相异号
            } else {
                self.calc_lnfp(self.Z) // 气相同号
            }
        }
    }
}
#[allow(non_snake_case)]
impl Flash for Pr {
    fn tp_flash(&mut self, T: f64, p: f64) -> Result<(), MyErr> {
        self.calc_AB(T, p);
        if !self.calc_root() {
            let lnfpg = self.calc_lnfp(self.Zg);
            let lnfpl = self.calc_lnfp(self.Zl);
            if lnfpg < lnfpl {
                self.Z = self.Zg;
            } else {
                self.Z = self.Zl;
            }
        }
        self.is_single_phase = true; // TPflash 得到单相区
        Ok(())
    }
    fn t_flash(&mut self, T: f64) -> Result<(), MyErr> {
        // 1. 寻找下边界
        let mut ps_min = 1.0; // 近零压力 1Pa
        self.calc_AB(T, ps_min);
        let mut lnfp_pmin = self.calc_diff_lnfpgl();
        // 2. 寻找上边界
        let mut ps_max = self.pc - 1.0; // 低于临界压力 1Pa
        self.calc_AB(T, ps_max);
        let mut lnfp_pmax = self.calc_diff_lnfpgl();
        // 3. 判断是否有解
        if lnfp_pmin * lnfp_pmax > 0.0 {
            if T < 0.95 * self.Tc {
                return Err(MyErr::new("T too low"));
            } else {
                return Err(MyErr::new("T too high"));
            }
        }
        // 4. 二分法计算
        // 随着压力的增加 lnfp 单调增加 不要问原因 因为我也不理解
        let mut counter = 0;
        let mut lnfp;
        let mut ps = (ps_min + ps_max) / 2.0;
        loop {
            if counter == 1000 {
                return Err(MyErr::new("divergence"));
            } else {
                counter += 1;
            }
            self.calc_AB(T, ps);
            lnfp = self.calc_diff_lnfpgl();
            if lnfp.abs() < 1E-6 {
                return Ok(());
            } else if lnfp * lnfp_pmin < 0.0 {
                ps_max = ps;
                lnfp_pmax = lnfp;
            } else if lnfp * lnfp_pmax < 0.0 {
                ps_min = ps;
                lnfp_pmin = lnfp;
            }
            ps = (ps_min + ps_max) / 2.0; // 二分饱和蒸汽压
        }
    }
}
impl Prop for Pr {
    fn T(&self) -> Result<f64, MyErr> {
        Ok(self.T)
    }
    fn p(&self) -> Result<f64, MyErr> {
        if self.is_single_phase {
            Ok(self.P)
        } else {
            Err(MyErr::new("no p in double phase"))
        }
    }
    fn rho(&self) -> Result<f64, MyErr> {
        if self.is_single_phase {
            Ok(self.P / self.Z / self.R / self.T)
        } else {
            Err(MyErr::new("no rho in double phase"))
        }
    }
    fn ps(&self) -> Result<f64, MyErr> {
        if self.is_single_phase {
            Err(MyErr::new("no ps in single phase"))
        } else {
            Ok(self.P)
        }
    }
    fn rhogs(&self) -> Result<f64, MyErr> {
        if self.is_single_phase {
            Err(MyErr::new("no rhogs in single phase"))
        } else {
            Ok(self.P / self.Zg / self.R / self.T)
        }
    }
    fn rhols(&self) -> Result<f64, MyErr> {
        if self.is_single_phase {
            Err(MyErr::new("no rhols in single phase"))
        } else {
            Ok(self.P / self.Zl / self.R / self.T)
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pr() {
        let mut n2 = Pr::new_fluid(126.19, 3.3958E6, 8.314, 0.0372);
        n2.tp_flash(300.0, 0.1E6).unwrap();
    }
}
