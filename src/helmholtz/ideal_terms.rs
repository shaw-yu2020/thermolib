use serde::{Deserialize, Serialize};
///
/// 理想气体亥姆霍兹方程
///
#[derive(Serialize, Deserialize, Debug)]
struct IdealPolynomialTerm {
    a: f64,
    t: f64,
}
#[derive(Serialize, Deserialize, Debug)]
struct IdealPlankEinsteinTerm {
    v: f64,
    b: f64,
}
