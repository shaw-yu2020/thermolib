#[allow(non_snake_case)]
pub fn rho_pr(T: f64, p: f64, Tc: f64, pc: f64, R: f64, omega: f64) -> f64 {
    let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega.powi(2);
    let a = 0.45724 * R.powi(2) * Tc.powi(2) / pc * (1.0 + kappa * (1.0 - (T / Tc).sqrt())).powi(2);
    let b = 0.07780 * R * Tc / pc;
    let A = a * p / R.powi(2) / T.powi(2);
    let B = b * p / R / T;
    let b = -1.0 + B;
    let c = A - 3.0 * B.powi(2) - 2.0 * B;
    let d = -A * B + B.powi(2) + B.powi(3);
    let A = b.powi(2) - 3.0 * c;
    let B = b * c - 9.0 * d;
    let C = c.powi(2) - 3.0 * b * d;
    let (Zv, Zl): (f64, f64);
    let Delta = B.powi(2) - 4.0 * A * C;
    if Delta.is_sign_negative() {
        let theta3 = ((2.0 * A * b - 3.0 * B) / (2.0 * A * A.sqrt())).acos() / 3.0;
        let x1 = (-b - 2.0 * A.sqrt() * theta3.cos()) / 3.0;
        let x2 = (-b + A.sqrt() * (theta3.cos() + 3_f64.sqrt() * theta3.sin())) / 3.0;
        let x3 = (-b + A.sqrt() * (theta3.cos() - 3_f64.sqrt() * theta3.sin())) / 3.0;
        Zv = x1.max(x2).max(x3);
        Zl = x1.min(x2).min(x3);
    } else {
        let Y1 = A * b + 1.5 * (-B + Delta.sqrt());
        let Y2 = A * b + 1.5 * (-B - Delta.sqrt());
        Zv = (-b - (Y1.cbrt() + Y2.cbrt())) / 3.0;
        Zl = -0.0;
    }
    if Zl.is_sign_negative() {
        p / (Zv * R * T)
    } else {
        let lnfpg = (Zv - 1.0)
            - (Zv - B).ln()
            - A / (2.0 * 2_f64.sqrt() * B) * ((Zv + 2.414 * B) / (Zv - 0.414 * B)).ln();
        let lnfpl = (Zl - 1.0)
            - (Zl - B).ln()
            - A / (2.0 * 2_f64.sqrt() * B) * ((Zl + 2.414 * B) / (Zl - 0.414 * B)).ln();
        if lnfpg < lnfpl {
            p / (Zv * R * T)
        } else {
            p / (Zl * R * T)
        }
    }
}
