// 偏导数
pub enum AlphaDtauDdeltaPlus {
    D00, // 温度0 密度0
    D01, // 温度0 密度1
    D02, // 温度0 密度2
    D10, // 温度1 密度0
    D20, // 温度2 密度0
}
// 辅助方程
mod ancillary_equations;
mod ideal_helmholtz_equation;
mod residual_helmholtz_equation;
// 亥姆霍兹方程
mod helmholtz_pure;
pub use helmholtz_pure::read_json;
pub use helmholtz_pure::HelmholtzPure;
