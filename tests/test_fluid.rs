/// Integration test
mod internals;
use internals::Value;
use thermolib::Helmholtz;
#[test]
#[allow(non_snake_case)]
fn test_SO2() {
    let mut SO2: Helmholtz = Helmholtz::read_json("SO2.json").expect("no SO2.json");
    let vec_value = vec![
        Value::new_short(250.0, 23.6e3, 12.2955804e6, 53.1514, 86.0982, 1130.24),
        Value::new_short(400.0, 16e3, 8.079379e6, 51.8705, 117.691, 449.618),
        Value::new_short(431.0, 8.078e3, 7.9343772e6, 64.5073, 19127.4, 168.147),
        Value::new_short(250.0, 0.0, 0.0, 29.8406, 38.1551, 203.682),
        Value::new_short(420.0, 1e3, 2.9365903e6, 40.8928, 59.5297, 234.103),
        Value::new_short(450.0, 11e3, 12.1084452e6, 54.787, 222.083, 250.095),
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
            14.88551e3,
            f64::INFINITY,
        ),
        Value::new_long(
            225.0,
            7.8e3,
            37.55722e6,
            192.9387,
            246.576,
            758.5512,
            -8.747458e3,
            -55.99954,
        ),
        Value::new_long(
            360.0, 5.2e3, 3.12811e6, 223.0894, 303.2828, 226.8389, 24.85197e3, 77.45943,
        ),
        Value::new_long(
            387.0, 2.637e3, 2.35539e6, 236.7771, 8976.589, 49.32548, 37.43821e3, 111.3057,
        ),
        Value::new_long(
            380.0,
            0.35e3,
            0.9312025e6,
            218.1084,
            236.858,
            99.66618,
            44.83037e3,
            135.8261,
        ),
        Value::new_long(
            400.0, 3.6e3, 3.513083e6, 233.3552, 437.7846, 89.44035, 38.4044e3, 112.8369,
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
            15.00556e3,
            f64::INFINITY,
        ),
        Value::new_long(
            250.0,
            6.5e3,
            45.74829e6,
            219.3241,
            268.4391,
            769.7973,
            -10.43608e3,
            -64.99865,
        ),
        Value::new_long(
            390.0, 4.2e3, 1.496384e6, 273.9917, 375.316, 182.6921, 29.10127e3, 83.3627,
        ),
        Value::new_long(
            421.5, 2.17e3, 2.083314e6, 302.6768, 15207.46, 41.76442, 44.91926e3, 121.6715,
        ),
        Value::new_long(
            410.0, 0.3e3, 0.841555e6, 273.0757, 294.6374, 91.41883, 51.90776e3, 143.3826,
        ),
        Value::new_long(
            450.0, 3e3, 4.15919e6, 297.3112, 431.0872, 99.09973, 51.36937e3, 134.6778,
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
            9.95556e3,
            f64::INFINITY,
        ),
        Value::new_long(
            260.0,
            5.5e3,
            28.03371e6,
            270.3971,
            329.6023,
            730.8597,
            -21.42844e3,
            -91.21752,
        ),
        Value::new_long(
            410.0,
            3.7e3,
            0.9573522e6,
            336.7461,
            435.6546,
            181.2565,
            31.64675e3,
            85.06658,
        ),
        Value::new_long(
            448.5, 1.825e3, 1.758863e6, 363.5137, 16301.72, 36.5901, 53.05034e3, 133.9994,
        ),
        Value::new_long(
            430.0,
            0.23e3,
            0.6728015e6,
            330.4599,
            352.0829,
            85.33926,
            58.08506e3,
            150.4846,
        ),
        Value::new_long(
            460.0, 2.7e3, 2.671262e6, 358.05, 541.4677, 81.00252, 53.96219e3, 135.1519,
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
            265.0,
            8.3e3,
            34.55032e6,
            137.1533,
            176.3509,
            1421.979,
            -10.80643e3,
            -50.64068,
        ),
        Value::new_long(
            465.0, 5.5e3, 3.513147e6, 208.5403, 278.1085, 394.104, 30.76257e3, 75.55693,
        ),
        Value::new_long(
            506.5, 2.78e3, 3.20769e6, 256.1725, 23439.18, 80.18092, 47.393e3, 109.6177,
        ),
        Value::new_long(
            265.0,
            0.0,
            0.0,
            118.0127,
            126.3272,
            165.4369,
            18.44498e3,
            f64::INFINITY,
        ),
        Value::new_long(
            485.0, 0.37e3, 1.215479e6, 207.5595, 226.9415, 181.491, 53.55585e3, 127.1686,
        ),
        Value::new_long(
            525.0, 3.8e3, 4.791053e6, 230.9643, 395.2492, 165.7838, 49.2901e3, 112.3864,
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
            250.0,
            8.3e3,
            36.16117e6,
            131.5927,
            168.4855,
            1433.881,
            -10.34456e3,
            -52.29832,
        ),
        Value::new_long(
            450.0, 5.5e3, 3.84858e6, 205.7891, 268.0357, 392.0591, 29.46773e3, 74.87199,
        ),
        Value::new_long(
            490.5, 2.78e3, 3.1618e6, 246.6004, 19167.87, 82.45625, 45.45018e3, 108.8356,
        ),
        Value::new_long(
            250.0,
            0.0,
            0.0,
            113.9448,
            122.2593,
            160.8752,
            16.76637e3,
            f64::INFINITY,
        ),
        Value::new_long(
            470.0, 0.37e3, 1.192083e6, 204.3634, 222.5305, 180.5417, 51.31643e3, 126.2985,
        ),
        Value::new_long(
            510.0, 3.8e3, 4.816117e6, 229.2782, 393.3932, 162.7817, 47.55079e3, 112.0508,
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
            260.0,
            8.3e3,
            33.95611e6,
            137.4007,
            171.8992,
            1402.545,
            -10.45592e3,
            -50.03608,
        ),
        Value::new_long(
            460.0, 5.5e3, 3.333684e6, 205.0061, 274.9322, 394.4134, 30.36433e3, 75.59279,
        ),
        Value::new_long(
            501.0, 2.8e3, 3.180216e6, 227.2418, 18442.03, 86.22365, 46.71702e3, 109.3744,
        ),
        Value::new_long(
            260.0,
            0.0,
            0.0,
            118.9821,
            127.2965,
            163.8248,
            17.84952e3,
            f64::INFINITY,
        ),
        Value::new_long(
            480.0, 0.37e3, 1.22339e6, 205.9278, 224.3395, 183.0561, 52.98636e3, 127.3124,
        ),
        Value::new_long(
            520.0, 3.8e3, 4.763216e6, 229.0952, 402.1148, 161.0, 48.68907e3, 112.3156,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut db);
    }
}
#[test]
#[allow(non_snake_case)]
fn test_NH3() {
    let mut NH3: Helmholtz = Helmholtz::read_json("NH3.json").expect("no NH3.json");
    let vec_value = vec![
        Value::new_long(
            200.0, 50e3, 555.5631e6, 56.15976, 74.22262, 2940.432, 9.541211e3, -11.61817,
        ),
        Value::new_long(
            250.0, 42e3, 108.7263e6, 50.25907, 72.04765, 2049.168, 5.715934e3, 14.26145,
        ),
        Value::new_long(
            320.0, 35e3, 27.63275e6, 46.95882, 78.92022, 1370.446, 9.831787e3, 36.0605,
        ),
        Value::new_long(
            380.0, 26e3, 8.174128e6, 47.82965, 125.1971, 668.5456, 15.45317e3, 53.81385,
        ),
        Value::new_long(
            405.0, 16.5e3, 11.25287e6, 77.13442, 4137.738, 247.2696, 20.13541e3, 65.30611,
        ),
        Value::new_long(
            // change pressure 11.45246e6 Pa to 11.452456
            406.0,
            13.696e3,
            11.452456e6,
            103.0508,
            52943.69,
            225.4458,
            21.3041e3,
            68.15535,
        ),
        Value::new_long(
            225.0,
            0.0,
            0.0,
            25.60265,
            33.91711,
            381.4709,
            26.38996e3,
            f64::INFINITY,
        ),
        Value::new_long(
            275.0,
            0.2e3,
            0.4270652e6,
            32.60471,
            45.39676,
            402.6253,
            27.46383e3,
            104.2956,
        ),
        Value::new_long(
            325.0, 0.95e3, 2.093496e6, 40.6866, 65.25439, 402.2104, 27.92681e3, 94.19375,
        ),
        Value::new_long(
            375.0, 3e3, 6.044096e6, 48.9189, 111.3738, 379.43, 27.45404e3, 86.36105,
        ),
        Value::new_long(
            405.0, 10e3, 11.2204e6, 77.15007, 2051.793, 286.2398, 23.20689e3, 72.89756,
        ),
        Value::new_long(
            520.0, 1.5e3, 5.96888e6, 36.52118, 49.38647, 538.8044, 36.23166e3, 106.6103,
        ),
        Value::new_long(
            520.0, 13.696e3, 35.83377e6, 43.66588, 98.36082, 555.4454, 28.71997e3, 80.49762,
        ),
        Value::new_long(
            520.0, 20e3, 56.3075e6, 43.81573, 87.32996, 741.1842, 26.48937e3, 73.89537,
        ),
        Value::new_long(
            620.0, 5.5e3, 24.38934e6, 41.29819, 61.81978, 588.0372, 38.40927e3, 100.0682,
        ),
        Value::new_long(
            620.0, 14e3, 58.52754e6, 43.50218, 75.46188, 692.5867, 34.49909e3, 87.78122,
        ),
        Value::new_long(
            620.0, 45e3, 1355.431e6, 50.78124, 64.36735, 3201.084, 50.78325e3, 54.81426,
        ),
    ];
    for value in vec_value.iter() {
        value.test(&mut NH3);
    }
}
