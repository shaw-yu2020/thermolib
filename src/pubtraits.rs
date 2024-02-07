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
                                       /*
                                       fn Z(&self) -> Result<f64, MyErr>; // 压缩因子
                                       fn cv(&self) -> Result<f64, MyErr>; // 定容比热
                                       fn cp(&self) -> Result<f64, MyErr>; // 定压比热
                                       fn w(&self) -> Result<f64, MyErr>; // 声速
                                       fn u(&self) -> Result<f64, MyErr>; // 比内能
                                       fn h(&self) -> Result<f64, MyErr>; // 比焓
                                       fn s(&self) -> Result<f64, MyErr>; // 比熵
                                       fn a(&self) -> Result<f64, MyErr>; // 比亥姆霍兹能
                                       fn g(&self) -> Result<f64, MyErr>; // 比吉布斯能
                                       */
    fn ps(&self) -> Result<f64, MyErr>; // 饱和蒸汽压
    fn rhogs(&self) -> Result<f64, MyErr>; // 饱和气相密度
    fn rhols(&self) -> Result<f64, MyErr>; // 饱和液相密度
}
