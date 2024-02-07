///
/// pubtraits
///
// 记录错误
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
    fn tp_flash(&mut self, T: f64, p: f64) -> Result<(), MyErr>; // 温度压力
    fn t_flash(&mut self, T: f64) -> Result<(), MyErr>; // 饱和温度
}
// 得到结果
#[allow(non_snake_case)]
pub trait Prop {
    fn T(&self) -> Result<f64, MyErr>; // 温度
    fn rho(&self) -> Result<f64, MyErr>; // 密度
    fn p(&self) -> Result<f64, MyErr>; // 压力
    fn ps(&self) -> Result<f64, MyErr>; // 饱和蒸汽压
    fn rhogs(&self) -> Result<f64, MyErr>; // 饱和气相密度
    fn rhols(&self) -> Result<f64, MyErr>; // 饱和液相密度
}
