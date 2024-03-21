///
/// 测试计算程序的可靠性
///
mod internals;
#[test]
fn test_so2() {
    let mut so2: thermolib::HelmholtzPure =
        thermolib::HelmholtzPure::read_json("SO2.json").expect("no SO2.json");
    let vec_vvs = vec![
        // 第一个数据点 来自GAO2016
        internals::new_value(250.0, 23600.0, 12295580.4, 53.1514, 86.0982, 1130.24),
        // 第二个数据点 来自GAO2016
        internals::new_value(400.0, 16000.0, 8079379.0, 51.8705, 117.691, 449.618),
        // 第三个数据点 来自GAO2016
        internals::new_value(431.0, 8078.0, 7934377.2, 64.5073, 19127.4, 168.147),
        // 第四个数据点 来自GAO2016
        internals::new_value(250.0, 0.0, 0.0, 29.8406, 38.1551, 203.682),
        // 第五个数据点 来自GAO2016
        internals::new_value(420.0, 1000.0, 2936590.3, 40.8928, 59.5297, 234.103),
        // 第六个数据点 来自GAO2016
        internals::new_value(450.0, 11000.0, 12108445.2, 54.787, 222.083, 250.095),
    ];
    internals::verify_fluid(&mut so2, &vec_vvs);
}
