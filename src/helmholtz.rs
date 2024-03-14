// 理想气体无量纲亥姆霍兹自由能 Alpha0
// 对无量纲温度的导数 Dtau
pub enum Alpha0Dtau {
    D0, // 零阶导数 实际是：tau^(0)*Alpha0Dtau0 -ln(delta)
    D1, // 一阶导数 实际是：tau^(1)*Alpha0Dtau1
    D2, // 二阶导数 实际是：tau^(2)*Alpha0Dtau2
}
// 无量纲亥姆霍兹自由能 Alpha
// 对无量纲温度 tau
// 和无量纲密度 delta
// 的偏导数
pub enum AlphaDD {
    D00, // 温度0 密度0
    D01, // 温度0 密度1
    D02, // 温度0 密度2
    D10, // 温度1 密度0
    D11, // 温度1 密度1
    D20, // 温度2 密度0
}
// 流体的热力学物性
pub enum ThermoProp {
    T,  // 温度
    D,  // 密度
    Z,  // 压缩因子
    P,  // 压力
    CV, // 定容比热
    CP, // 定压比热
    W,  // 声速
    S,  // 比熵
    U,  // 比内能
    H,  // 比焓
    A,  // 比亥姆霍兹能
    G,  // 比吉布斯能
}
// 辅助方程
mod ancillary_equations;
mod ideal_helmholtz_equation;
mod residual_helmholtz_equation;
// 实际气体亥姆霍兹方程
mod real_helmholtz_equation;
// 亥姆霍兹方程
mod helmholtz_pure;
pub use helmholtz_pure::read_json;
pub use helmholtz_pure::HelmholtzPure;
