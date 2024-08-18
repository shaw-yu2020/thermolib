/// Reference:
/// Richard Brent,
/// Algorithms for Minimization Without Derivatives.
pub fn brent_zero<F>(mut f: F, a: f64, b: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    let macheps = f64::EPSILON.sqrt();
    let t = 10.0 * f64::EPSILON.sqrt();
    let (mut a, mut fa) = (a, f(a));
    let (mut b, mut fb) = (b, f(b));
    if (fa * fb) > 0.0 {
        return f64::NAN;
    }
    let (mut c, mut fc) = (a, fa);
    let (mut d, mut e) = (b - a, b - a);
    let (mut tol, mut m, mut p, mut q, mut r, mut s);
    loop {
        if fc.abs() < fb.abs() {
            (a, fa) = (b, fb);
            (b, fb) = (c, fc);
            (c, fc) = (a, fa);
        }
        tol = 2.0 * macheps * b.abs() + t;
        m = 0.5 * (c - b);
        if m.abs() <= tol || fb == 0.0 {
            return b;
        }
        if e.abs() < tol || fa.abs() <= fb.abs() {
            e = m;
            d = e;
        } else {
            s = fb / fa;
            if a == c {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if 0.0 < p {
                q = -q;
            } else {
                p = -p;
            }
            s = e;
            e = d;
            if 2.0 * p < 3.0 * m * q - (tol * q).abs() && p < (0.5 * s * q).abs() {
                d = p / q;
            } else {
                e = m;
                d = e;
            }
        }
        (a, fa) = (b, fb);
        b += if tol < d.abs() {
            d
        } else if 0.0 < m {
            tol
        } else {
            -tol
        };
        fb = f(b);
        if (0.0 < fb && 0.0 < fc) || (fb.is_sign_negative() && fc.is_sign_negative()) {
            (c, fc) = (a, fa);
            e = b - a;
            d = e;
        }
    }
}
/// unit test
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_algorithm() {
        let x0 = brent_zero(|x: f64| x.sin() - 0.5 * x, 1.0, 2.0);
        let (x1, digits) = (1.89549, 10_f64.powi(5));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| 2.0 * x - (-x).exp(), 0.0, 1.0);
        let (x1, digits) = (0.351734, 10_f64.powi(6));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| x * (-x).exp(), -1.0, 0.5);
        let (x1, digits) = (-4.03522e-10, 10_f64.powi(15));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| x.exp() - 1.0 / 100.0 / x / x, 0.0001, 20.0);
        let (x1, digits) = (0.0953446, 10_f64.powi(7));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| (x + 3.0) * (x - 1.0) * (x - 1.0), -5.0, 20.0);
        let (x1, digits) = (-3.0, 10_f64.powi(0));
        assert_eq!(x1, (x0 * digits).round() / digits);
    }
}
