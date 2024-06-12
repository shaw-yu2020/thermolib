/// Integration test
use serde::Deserialize;
use thermolib::Helmholtz;
#[test]
#[allow(non_snake_case)]
fn test_fluids() {
    let entries = std::fs::read_dir(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("res"))
        .expect("Failed to read directory");
    let mut vec_fluid_json: Vec<String> = Vec::new();
    for entry in entries.flatten() {
        vec_fluid_json.push(entry.file_name().to_string_lossy().to_string())
    }
    for fluid_json in vec_fluid_json.iter() {
        let mut file_json = std::fs::File::open(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("res")
                .join(fluid_json),
        )
        .unwrap();
        let mut str_json = String::new();
        let _ = std::io::Read::read_to_string(&mut file_json, &mut str_json);
        let mut fluid: Helmholtz = serde_json::from_str(&str_json).unwrap();
        let val: serde_json::Value = serde_json::from_str(&str_json).unwrap();
        let Tt: f64 = val["Tt"].as_f64().unwrap();
        let Tc: f64 = val["Tc"].as_f64().unwrap();
        let mut T = Tt.floor();
        loop {
            if let Err(_) = fluid.t_flash(T) {
                if Tc.fract() == 0.0 {
                    assert_eq!(Tc + 1.0, T);
                } else {
                    assert_eq!(Tc.ceil(), T);
                }
                break;
            }
            T += 1.0;
        }
        let values: ValueS = serde_json::from_str(&str_json).unwrap();
        for value in values.values.iter() {
            value.test(&mut fluid);
        }
    }
}
#[derive(Deserialize)]
pub struct ValueS {
    #[serde(default)]
    pub values: Vec<Value>,
}
#[derive(Deserialize)]
#[allow(non_snake_case)]
pub struct Value {
    T: f64,
    rho: f64,
    p: f64,
    cv: f64,
    cp: f64,
    w: f64,
    #[serde(default)]
    h: f64,
    #[serde(default)]
    s: f64,
}
impl Value {
    pub fn test(&self, fluid: &mut Helmholtz) {
        // test td_flash
        fluid.td_flash(self.T, self.rho).unwrap();
        compare_eq(self.p, fluid.p());
        compare_eq(self.cv, fluid.cv().unwrap());
        compare_eq(self.cp, fluid.cp().unwrap());
        compare_eq(self.w, fluid.w().unwrap());
        if 0.0 != self.h {
            compare_eq(self.h, fluid.h().unwrap());
        }
        if 0.0 != self.s {
            if 1e16 == self.s {
                compare_eq(f64::INFINITY, fluid.s().unwrap());
            } else {
                compare_eq(self.s, fluid.s().unwrap());
            }
        }
        // test tp_flash
        if self.p != 0.0 {
            fluid.tp_flash(self.T, self.p).unwrap();
            compare_eq(self.rho, fluid.rho().unwrap());
        }
    }
}
fn compare_eq(f64_short: f64, f64_long: f64) {
    if (f64_short / f64_long - f64_long / f64_short).abs() / 2.0 * 1e12 < 1.0 {
        return;
    } // check relative deviation.
    let string_short = f64_short.to_string();
    let mut length = match string_short.rfind(".") {
        Some(index) => (string_short.len() - index - 1) as i32,
        None => 0,
    };
    if 0 == length && 0.0 != f64_short && 0.0 != f64_long {
        while 0.0 == f64_short % 10_f64.powi(-length) {
            length -= 1;
        }
        length += 1;
    };
    let round_long = (f64_long * 10_f64.powi(length)).round() / 10_f64.powi(length);
    assert_eq!(f64_short, round_long);
}
