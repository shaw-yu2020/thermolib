/// Reference:
/// Richard Brent,
/// Algorithms for Minimization Without Derivatives.
pub fn brent_zero<F>(mut f: F, a: f64, b: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    let macheps = f64::EPSILON.sqrt();
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
        tol = macheps * (2.0 * b.abs() + 10.0);
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
/// Romberg Numerical Differentiation
/// Reference:
/// C. J. F. RIDDERS
/// Accurate computation of F'(x) and F'(x) F"(x)
/// Adv. Eng. Software, 1982, Vol. 4, No. 2, 75-76.
const N_DIM: usize = 11;
pub fn romberg_diff<F>(mut f: F, x: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    let tol = f64::EPSILON.sqrt() * 2.0;
    let mut fx: [[f64; N_DIM]; N_DIM] = [[0.0; N_DIM]; N_DIM];
    let mut h: f64 = if x > 1.0 { 0.01 * x } else { 0.1 * x };
    fx[0][0] = (f(x + h) - f(x - h)) / 2.0 / h;
    for i in 1..N_DIM {
        h /= 2.0;
        fx[0][i] = (f(x + h) - f(x - h)) / 2.0 / h;
        for j in 1..i + 1 {
            fx[j][i] = (fx[j - 1][i] * 4_f64.powi(j as i32) - fx[j - 1][i - 1])
                / (4_f64.powi(j as i32) - 1.0);
        }
        if (fx[i][i] / fx[i - 1][i - 1] - fx[i - 1][i - 1] / fx[i][i]).abs() < tol {
            return fx[i][i];
        }
    }
    fx[N_DIM - 1][N_DIM - 1]
}
/// unit test
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_algorithms() {
        // test brent_zero
        let x0 = brent_zero(|x: f64| x.sin() - x / 2.0, 1.0, 2.0);
        let (x1, digits) = (1.8954942796e+00, 10_f64.powi(10));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| 2.0 * x - (-x).exp(), 0.0, 1.0);
        let (x1, digits) = (3.5173365631e-01, 10_f64.powi(11));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| x * (-x).exp(), -1.0, 0.5);
        let (x1, digits) = (-4.0352160429e-10, 10_f64.powi(20));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| x.exp() - 1.0 / 100.0 / x / x, 0.0001, 20.0);
        let (x1, digits) = (9.5344620258e-02, 10_f64.powi(12));
        assert_eq!(x1, (x0 * digits).round() / digits);
        let x0 = brent_zero(|x: f64| (x + 3.0) * (x - 1.0) * (x - 1.0), -5.0, 2.0);
        let (x1, digits) = (-3.0000000000e+00, 10_f64.powi(10));
        assert_eq!(x1, (x0 * digits).round() / digits);
        // test romberg_diff
        let (df0, digits) = (140.7377356, 10_f64.powi(7)); // real = 140.7377355
        let df1 = romberg_diff(|x: f64| x.exp() / (x.sin() - x.powi(2)), 1.0);
        assert_eq!(df0, (df1 * digits).round() / digits);
        let (df0, digits) = (2.718281828, 10_f64.powi(9));
        let df1 = romberg_diff(|x: f64| x.exp(), 1.0);
        assert_eq!(df0, (df1 * digits).round() / digits);
        let (df0, digits) = (22026.46579, 10_f64.powi(5));
        let df1 = romberg_diff(|x: f64| x.exp(), 10.0);
        assert_eq!(df0, (df1 * digits).round() / digits);
    }
}
