///
/// 亥姆霍兹方程
///
#[allow(non_snake_case)]
struct CriticalPoint {
    Tc: f64,   // 临界温度
    rhoc: f64, // 临界密度
    pc: f64,   // 临界压力
    M: f64,    // 摩尔质量
    R: f64,    // 偏心因子
}
struct ResidualPolynomialTerm {
    n: f64,
    d: f64,
    t: f64,
}
struct ResidualExponentialTerm {
    n: f64,
    d: f64,
    t: f64,
    l: f64,
}
struct ResidualGaussBellTerm {
    n: f64,
    d: f64,
    t: f64,
    eta: f64,
    epsilon: f64,
    beta: f64,
    gamma: f64,
}
