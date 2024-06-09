/// Internal module
use thermolib::Helmholtz;
fn compare_eq(f64_short: f64, f64_long: f64) {
    if (f64_short / f64_long - f64_long / f64_short).abs() / 2.0 * 1e12 < 1.0 {
        return;
    } // check relative deviation.
    let string_short = f64_short.to_string();
    let length = match string_short.rfind(".") {
        Some(index) => (string_short.len() - index - 1) as i32,
        None => 0,
    };
    let round_long = (f64_long * 10_f64.powi(length)).round() / 10_f64.powi(length);
    assert_eq!(f64_short, round_long);
}
#[allow(non_snake_case)]
pub struct Value {
    T: f64,
    rho: f64,
    p: f64,
    cv: f64,
    cp: f64,
    w: f64,
    h: Option<f64>,
    s: Option<f64>,
}
#[allow(non_snake_case)]
impl Value {
    pub fn new_short(T: f64, rho: f64, p: f64, cv: f64, cp: f64, w: f64) -> Value {
        Value {
            T,
            rho,
            p,
            cv,
            cp,
            w,
            h: None,
            s: None,
        }
    }
    pub fn new_long(T: f64, rho: f64, p: f64, cv: f64, cp: f64, w: f64, h: f64, s: f64) -> Value {
        Value {
            T,
            rho,
            p,
            cv,
            cp,
            w,
            h: Some(h),
            s: Some(s),
        }
    }
    pub fn test(&self, fluid: &mut Helmholtz) {
        // test td_flash
        fluid.td_flash(self.T, self.rho).unwrap();
        compare_eq(self.p, fluid.p());
        compare_eq(self.cv, fluid.cv().unwrap());
        compare_eq(self.cp, fluid.cp().unwrap());
        compare_eq(self.w, fluid.w().unwrap());
        if let Some(h) = self.h {
            compare_eq(h, fluid.h().unwrap());
        }
        if let Some(s) = self.s {
            compare_eq(s, fluid.s().unwrap());
        }
        // test tp_flash
        if self.p != 0.0 {
            fluid.tp_flash(self.T, self.p).unwrap();
            compare_eq(self.rho, fluid.rho().unwrap());
        }
    }
}
