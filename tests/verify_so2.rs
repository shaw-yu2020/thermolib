///
/// 测试计算程序的可靠性
///
use thermolib::HelmholtzPure;
use thermolib::ThermoProp;
#[test]
fn test_so2() {
    let so2: HelmholtzPure = thermolib::read_json("SO2.json");
    let vec_vvs = vec![
        VerificationValue {
            T: 250.0,
            D: 23600.0,
            P: 12295580.4,
            CV: 53.1514,
            CP: 86.0982,
            W: 1130.24,
        }, // 第一个数据点 来自GAO2016
        VerificationValue {
            T: 400.0,
            D: 16000.0,
            P: 8079379.0,
            CV: 51.8705,
            CP: 117.691,
            W: 449.618,
        }, // 第二个数据点 来自GAO2016
        VerificationValue {
            T: 431.0,
            D: 8078.0,
            P: 7934377.2,
            CV: 64.5073,
            CP: 19127.4,
            W: 168.147,
        }, // 第三个数据点 来自GAO2016
        VerificationValue {
            T: 250.0,
            D: 0.0,
            P: 0.0,
            CV: 29.8406,
            CP: 38.1551,
            W: 203.682,
        }, // 第四个数据点 来自GAO2016
        VerificationValue {
            T: 420.0,
            D: 1000.0,
            P: 2936590.3,
            CV: 40.8928,
            CP: 59.5297,
            W: 234.103,
        }, // 第五个数据点 来自GAO2016
        VerificationValue {
            T: 450.0,
            D: 11000.0,
            P: 12108445.2,
            CV: 54.787,
            CP: 222.083,
            W: 250.095,
        }, // 第六个数据点 来自GAO2016
    ];
    verify(&so2, &vec_vvs);
}
#[allow(non_snake_case)]
struct VerificationValue {
    T: f64,
    D: f64,
    P: f64,
    CV: f64,
    CP: f64,
    W: f64,
}
fn verify(fluid: &HelmholtzPure, values: &Vec<VerificationValue>) {
    for value in values.iter() {
        compare_eq(value.P, fluid.calc(ThermoProp::P, value.T, value.D));
        compare_eq(value.CV, fluid.calc(ThermoProp::CV, value.T, value.D));
        compare_eq(value.CP, fluid.calc(ThermoProp::CP, value.T, value.D));
        compare_eq(value.W, fluid.calc(ThermoProp::W, value.T, value.D));
    }
}
fn compare_eq(f64_short: f64, f64_long: f64) {
    let string_short = f64_short.to_string();
    let length = match string_short.rfind(".") {
        Some(index) => (string_short.len() - index - 1) as i32,
        None => 0,
    };
    let round_long = (f64_long * 10_f64.powi(length)).round() / 10_f64.powi(length);
    assert_eq!(f64_short, round_long);
}
