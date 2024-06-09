/// Integration test
mod internals;
use internals::Value;
use thermolib::Helmholtz;
#[test]
#[allow(non_snake_case)]
fn test_SO2() {
    let mut SO2: Helmholtz = Helmholtz::read_json("SO2.json").expect("no SO2.json");
    let vec_value = vec![
        Value::new_short(250.0, 23600.0, 12295580.4, 53.1514, 86.0982, 1130.24),
        Value::new_short(400.0, 16000.0, 8079379.0, 51.8705, 117.691, 449.618),
        Value::new_short(431.0, 8078.0, 7934377.2, 64.5073, 19127.4, 168.147),
        Value::new_short(250.0, 0.0, 0.0, 29.8406, 38.1551, 203.682),
        Value::new_short(420.0, 1000.0, 2936590.3, 40.8928, 59.5297, 234.103),
        Value::new_short(450.0, 11000.0, 12108445.2, 54.787, 222.083, 250.095),
    ];
    for value in vec_value.iter() {
        value.test(&mut SO2);
    }
}
#[test]
#[allow(non_snake_case)]
fn test_C4F10() {
    let mut C4F10: Helmholtz = Helmholtz::read_json("C4F10.json").expect("no C4F10.json");
    let vec_value = vec![
        Value::new_long(
            225.0,
            0.0,
            0.0,
            171.226,
            179.5405,
            90.78029,
            14885.51,
            f64::INFINITY,
        ),
        Value::new_long(
            // change pressure 37557220.0 Pa to 37557221.0 Pa
            225.0, 7800.0, 37557221.0, 192.9387, 246.576, 758.5512, -8747.458, -55.99954,
        ),
        Value::new_long(
            360.0, 5200.0, 3128110.0, 223.0894, 303.2828, 226.8389, 24851.97, 77.45943,
        ),
        Value::new_long(
            387.0, 2637.0, 2355390.0, 236.7771, 8976.589, 49.32548, 37438.21, 111.3057,
        ),
        Value::new_long(
            380.0, 350.0, 931202.5, 218.1084, 236.858, 99.66618, 44830.37, 135.8261,
        ),
        Value::new_long(
            400.0, 3600.0, 3513083.0, 233.3552, 437.7846, 89.44035, 38404.44, 112.8369,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut C4F10);
    }
}
#[test]
#[allow(non_snake_case)]
fn test_C5F12() {
    let mut C5F12: Helmholtz = Helmholtz::read_json("C5F12.json").expect("no C5F12.json");
    let vec_value = vec![
        Value::new_long(
            250.0,
            0.0,
            0.0,
            199.2576,
            207.5721,
            86.70462,
            15005.56,
            f64::INFINITY,
        ),
        Value::new_long(
            // change pressure 45748290.0 Pa to 45748288.0 Pa
            250.0, 6500.0, 45748288.0, 219.3241, 268.4391, 769.7973, -10436.08, -64.99865,
        ),
        Value::new_long(
            390.0, 4200.0, 1496384.0, 273.9917, 375.316, 182.6921, 29101.27, 83.3627,
        ),
        Value::new_long(
            421.5, 2170.0, 2083314.0, 302.6768, 15207.46, 41.76442, 44919.26, 121.6715,
        ),
        Value::new_long(
            410.0, 300.0, 841555.0, 273.0757, 294.6374, 91.41883, 51907.76, 143.3826,
        ),
        Value::new_long(
            450.0, 3000.0, 4159190.0, 297.3112, 431.0872, 99.09973, 51369.37, 134.6778,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut C5F12);
    }
}
#[test]
#[allow(non_snake_case)]
fn test_C6F14() {
    let mut C6F14: Helmholtz = Helmholtz::read_json("C6F14.json").expect("no C6F14.json");
    let vec_value = vec![
        Value::new_long(
            260.0,
            0.0,
            0.0,
            244.6528,
            252.9673,
            81.3159,
            9955.56,
            f64::INFINITY,
        ),
        Value::new_long(
            // change pressure 28033710.0 Pa to 28033711.0 Pa
            260.0, 5500.0, 28033711.0, 270.3971, 329.6023, 730.8597, -21428.44, -91.21752,
        ),
        Value::new_long(
            410.0, 3700.0, 957352.2, 336.7461, 435.6546, 181.2565, 31646.75, 85.06658,
        ),
        Value::new_long(
            448.5, 1825.0, 1758863.0, 363.5137, 16301.72, 36.5901, 53050.34, 133.9994,
        ),
        Value::new_long(
            430.0, 230.0, 672801.5, 330.4599, 352.0829, 85.33926, 58085.06, 150.4846,
        ),
        Value::new_long(
            460.0, 2700.0, 2671262.0, 358.05, 541.4677, 81.00252, 53962.19, 135.1519,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut C6F14);
    }
}
#[test]
#[allow(non_snake_case)]
fn test_3MethylPentane() {
    let mut mp: Helmholtz =
        Helmholtz::read_json("3MethylPentane.json").expect("no 3MethylPentane.json");
    let vec_value = vec![
        Value::new_long(
            // change pressure 34550320.0 Pa to 34550316.0 Pa
            265.0, 8300.0, 34550316.0, 137.1533, 176.3509, 1421.979, -10806.43, -50.64068,
        ),
        Value::new_long(
            465.0, 5500.0, 3513147.0, 208.5403, 278.1085, 394.104, 30762.57, 75.55693,
        ),
        Value::new_long(
            506.5, 2780.0, 3207690.0, 256.1725, 23439.18, 80.18092, 47393.00, 109.6177,
        ),
        Value::new_long(
            265.0,
            0.0,
            0.0,
            118.0127,
            126.3272,
            165.4369,
            18444.98,
            f64::INFINITY,
        ),
        Value::new_long(
            485.0, 370.0, 1215479.0, 207.5595, 226.9415, 181.491, 53555.85, 127.1686,
        ),
        Value::new_long(
            525.0, 3800.0, 4791053.0, 230.9643, 395.2492, 165.7838, 49290.1, 112.3864,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut mp);
    }
}
#[test]
#[allow(non_snake_case)]
fn test_22DimethylButane() {
    let mut db: Helmholtz =
        Helmholtz::read_json("22DimethylButane.json").expect("no 22DimethylButane.json");
    let vec_value = vec![
        Value::new_long(
            // change pressure 36161170.0 Pa to 36161171.0 Pa
            250.0, 8300.0, 36161171.0, 131.5927, 168.4855, 1433.881, -10344.56, -52.29832,
        ),
        Value::new_long(
            450.0, 5500.0, 3848580.0, 205.7891, 268.0357, 392.0591, 29467.73, 74.87199,
        ),
        Value::new_long(
            490.5, 2780.0, 3161800.0, 246.6004, 19167.87, 82.45625, 45450.18, 108.8356,
        ),
        Value::new_long(
            250.0,
            0.0,
            0.0,
            113.9448,
            122.2593,
            160.8752,
            16766.37,
            f64::INFINITY,
        ),
        Value::new_long(
            470.0, 370.0, 1192083.0, 204.3634, 222.5305, 180.5417, 51316.43, 126.2985,
        ),
        Value::new_long(
            510.0, 3800.0, 4816117.0, 229.2782, 393.3932, 162.7817, 47550.79, 112.0508,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut db);
    }
}
#[test]
#[allow(non_snake_case)]
fn test_23DimethylButane() {
    let mut db: Helmholtz =
        Helmholtz::read_json("23DimethylButane.json").expect("no 23DimethylButane.json");
    let vec_value = vec![
        Value::new_long(
            // change pressure 33956110.0 Pa to 33956109.0 Pa
            260.0, 8300.0, 33956109.0, 137.4007, 171.8992, 1402.545, -10455.92, -50.03608,
        ),
        Value::new_long(
            460.0, 5500.0, 3333684.0, 205.0061, 274.9322, 394.4134, 30364.33, 75.59279,
        ),
        Value::new_long(
            501.0, 2800.0, 3180216.0, 227.2418, 18442.03, 86.22365, 46717.02, 109.3744,
        ),
        Value::new_long(
            260.0,
            0.0,
            0.0,
            118.9821,
            127.2965,
            163.8248,
            17849.52,
            f64::INFINITY,
        ),
        Value::new_long(
            480.0, 370.0, 1223390.0, 205.9278, 224.3395, 183.0561, 52986.36, 127.3124,
        ),
        Value::new_long(
            520.0, 3800.0, 4763216.0, 229.0952, 402.1148, 161.0, 48689.07, 112.3156,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut db);
    }
}
