///
/// pubtraits
///
// 记录错误
use std::fmt::Debug;
#[derive(Debug)]
pub struct MyErr {
    err: String,
}
impl MyErr {
    pub fn new(my_str: &str) -> MyErr {
        MyErr {
            err: String::from(my_str),
        }
    }
}
impl std::fmt::Display for MyErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.err)
    }
}
// 迭代计算
#[allow(non_snake_case)]
pub trait Flash {
    fn t_flash(&mut self, _T: f64) -> Result<(), MyErr>
    // 饱和温度
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn td_flash(&mut self, _T: f64, _rho: f64) -> Result<(), MyErr>
    // 温度密度
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn tp_flash(&mut self, _T: f64, _p: f64) -> Result<(), MyErr>
    // 温度压力
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
}
// 得到结果
#[allow(non_snake_case)]
pub trait Prop {
    fn T(&self) -> Result<f64, MyErr>
    // 温度
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn rho(&self) -> Result<f64, MyErr>
    // 密度
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn p(&self) -> Result<f64, MyErr>
    // 压力
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn cv(&self) -> Result<f64, MyErr>
    // 定容比热
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn cp(&self) -> Result<f64, MyErr>
    // 定压比热
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn w(&self) -> Result<f64, MyErr>
    // 声速
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn ps(&self) -> Result<f64, MyErr>
    // 饱和蒸汽压
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn rhogs(&self) -> Result<f64, MyErr>
    // 饱和气相密度
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
    fn rhols(&self) -> Result<f64, MyErr>
    // 饱和液相密度
    where
        Self: Debug,
    {
        Err(MyErr::new(&format!("no implementation for {:?}", self)))
    }
}
