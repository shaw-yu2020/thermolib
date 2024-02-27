use serde::{Deserialize, Serialize};
///
/// 剩余气体亥姆霍兹方程
///
#[derive(Serialize, Deserialize, Debug)]
struct ResidualPolynomialTerm {
    n: f64,
    d: f64,
    t: f64,
}
#[derive(Serialize, Deserialize, Debug)]
struct ResidualExponentialTerm {
    n: f64,
    d: f64,
    t: f64,
    l: f64,
}
#[derive(Serialize, Deserialize, Debug)]
struct ResidualGaussianTerm {
    n: f64,
    d: f64,
    t: f64,
    eta: f64,
    epsilon: f64,
    beta: f64,
    gamma: f64,
}
