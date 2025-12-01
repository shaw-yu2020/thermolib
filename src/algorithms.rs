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
/// y = c0 + c1 * x
/// s = min{ 0.5 * sum[ ( c0 + c1 * xi - yi )^2 ] }
///
/// ds/dc0 = sum[ ( c0 + c1 * xi - yi ) ]
///        = c0 * n + c1 * sum(xi) - sum(yi)
/// ds/dc1 = sum[ ( c0 + c1 * xi - yi ) * xi ]
///        = c0 * sum(xi) + c1 * sum(xi*xi) - sum(yi*xi)
///
/// | n       sum(xi)    | |c0| = | sum(yi)    |
/// | sum(xi) sum(xi*xi) | |c1| = | sum(yi*xi) |
///
/// ci = | n       sum(xi)    |
///      | sum(xi) sum(xi*xi) | = n*sum(xi*xi) - sum(xi)*sum(xi)
///
/// c0 = | sum(yi)    sum(xi)    |
///      | sum(yi*xi) sum(xi*xi) | / ci
/// c1 = | n       sum(yi)    |
///      | sum(xi) sum(yi*xi) | / ci
///
/// c0 = frac{ sum(yi)*sum(xi*xi) - sum(yi*xi)*sum(xi) }
///          { n*sum(xi*xi) - sum(xi)*sum(xi) }
/// c1 = frac{ n*sum(yi*xi) - sum(xi)*sum(yi) }
///          { n*sum(xi*xi) - sum(xi)*sum(xi) }
///
pub fn poly1fit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let (mut sum_xi, mut sum_xixi) = (0.0, 0.0);
    let (mut sum_yi, mut sum_yixi) = (0.0, 0.0);
    let num: f64 = (0..usize::min(x.len(), y.len()))
        .map(|i| {
            sum_xi += x[i];
            sum_yi += y[i];
            sum_xixi += x[i] * x[i];
            sum_yixi += y[i] * x[i];
        })
        .count() as f64;
    let ci = num * sum_xixi - sum_xi * sum_xi;
    if ci.abs() > 0.01 {
        (
            (sum_yi * sum_xixi - sum_yixi * sum_xi) / ci,
            (num * sum_yixi - sum_xi * sum_yi) / ci,
        )
    } else {
        let ata = vec![vec![num, sum_xi], vec![sum_xi, sum_xixi]];
        let atb = vec![sum_yi, sum_yixi];
        let result = solve_equ(ata, atb);
        (result[0], result[1])
    }
}
/// y = c0 + c1 * x + c2 * x * x
/// s = min{ 0.5 * sum[ ( c0 + c1 * xi + c2 * xi * xi - yi )^2 ] }
///
/// ds/dc0 = sum[ ( c0 + c1 * xi + c2 * xi * xi - yi ) ]
///        = c0 * n + c1 * sum(xi) + c2 * sum(xi*xi) - sum(yi)
/// ds/dc1 = sum[ ( c0 + c1 * xi + c2 * xi * xi - yi ) * xi ]
///        = c0 * sum(xi) + c1 * sum(xi*xi) + c2 * sum(xi*xi*xi) - sum(yi*xi)
/// ds/dc2 = sum[ ( c0 + c1 * xi + c2 * xi * xi - yi ) * xi * xi ]
///        = c0 * sum(xi*xi) + c1 * sum(xi*xi*xi) + c2 * sum(xi*xi*xi*xi) - sum(yi*xi*xi)
///
/// | n          sum(xi)       sum(xi*xi)       | |c0| = | sum(yi)       |
/// | sum(xi)    sum(xi*xi)    sum(xi*xi*xi)    | |c1| = | sum(yi*xi)    |
/// | sum(xi*xi) sum(xi*xi*xi) sum(xi*xi*xi*xi) | |c2| = | sum(yi*xi*xi) |
///
/// ci = | n          sum(xi)       sum(xi*xi)       |
///      | sum(xi)    sum(xi*xi)    sum(xi*xi*xi)    |
///      | sum(xi*xi) sum(xi*xi*xi) sum(xi*xi*xi*xi) |
///    = n*[ sum(xi*xi)*sum(xi*xi*xi*xi) - sum(xi*xi*xi)*sum(xi*xi*xi) ]
///    + sum(xi)*[ sum(xi*xi*xi)*sum(xi*xi) - sum(xi)*sum(xi*xi*xi*xi) ]
///    + sum(xi*xi)*[ sum(xi)*sum(xi*xi*xi) - sum(xi*xi)*sum(xi*xi) ]
///
/// c0 = | sum(yi)       sum(xi)       sum(xi*xi)       |
///      | sum(yi*xi)    sum(xi*xi)    sum(xi*xi*xi)    |
///      | sum(yi*xi*xi) sum(xi*xi*xi) sum(xi*xi*xi*xi) |
/// c1 = | n          sum(yi)       sum(xi*xi)       |
///      | sum(xi)    sum(yi*xi)    sum(xi*xi*xi)    |
///      | sum(xi*xi) sum(yi*xi*xi) sum(xi*xi*xi*xi) |
/// c2 = | n          sum(xi)       sum(yi)       |
///      | sum(xi)    sum(xi*xi)    sum(yi*xi)    |
///      | sum(xi*xi) sum(xi*xi*xi) sum(yi*xi*xi) |
///
pub fn poly2fit(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let (mut sum_x1, mut sum_x2, mut sum_x3, mut sum_x4) = (0.0, 0.0, 0.0, 0.0);
    let (mut sum_y1x0, mut sum_y1x1, mut sum_y1x2) = (0.0, 0.0, 0.0);
    let num: f64 = (0..usize::min(x.len(), y.len()))
        .map(|i| {
            let (x1, x2) = (x[i], x[i] * x[i]);
            sum_x1 += x1;
            sum_x2 += x2;
            sum_x3 += x1 * x2;
            sum_x4 += x2 * x2;
            sum_y1x0 += y[i];
            sum_y1x1 += y[i] * x1;
            sum_y1x2 += y[i] * x2;
        })
        .count() as f64;
    let ci = num * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
        + sum_x1 * (sum_x3 * sum_x2 - sum_x1 * sum_x4)
        + sum_x2 * (sum_x1 * sum_x3 - sum_x2 * sum_x2);
    if ci.abs() > 0.01 {
        (
            (sum_y1x0 * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
                + sum_y1x1 * (sum_x3 * sum_x2 - sum_x1 * sum_x4)
                + sum_y1x2 * (sum_x1 * sum_x3 - sum_x2 * sum_x2))
                / ci,
            (num * (sum_y1x1 * sum_x4 - sum_y1x2 * sum_x3)
                + sum_x1 * (sum_y1x2 * sum_x2 - sum_y1x0 * sum_x4)
                + sum_x2 * (sum_y1x0 * sum_x3 - sum_y1x1 * sum_x2))
                / ci,
            (num * (sum_x2 * sum_y1x2 - sum_x3 * sum_y1x1)
                + sum_x1 * (sum_x3 * sum_y1x0 - sum_x1 * sum_y1x2)
                + sum_x2 * (sum_x1 * sum_y1x1 - sum_x2 * sum_y1x0))
                / ci,
        )
    } else {
        let ata = vec![
            vec![num, sum_x1, sum_x2],
            vec![sum_x1, sum_x2, sum_x3],
            vec![sum_x2, sum_x3, sum_x4],
        ];
        let atb = vec![sum_y1x0, sum_y1x1, sum_y1x2];
        let result = solve_equ(ata, atb);
        (result[0], result[1], result[2])
    }
}
/// Solver -> Ax=b
fn solve_equ(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = a.len();
    let (mut index, mut scale);
    for k in 0..n {
        index = k;
        for i in k + 1..n {
            if a[i][k] > a[index][k] {
                index = i;
            }
        }
        if index > k {
            a.swap(k, index);
            b.swap(k, index);
        }
        for i in k + 1..n {
            scale = a[i][k] / a[k][k];
            for j in k..n {
                a[i][j] -= scale * a[k][j];
            }
            b[i] -= scale * b[k];
        }
    }
    let mut x = vec![0.0; b.len()];
    for k in (0..n).rev() {
        x[k] = (b[k]
            - a[k][k + 1..]
                .iter()
                .zip(x[k + 1..].iter())
                .map(|(aj, xj)| aj * xj)
                .sum::<f64>())
            / a[k][k];
    }
    x
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
        // test poly1fit -> y = 1 + 3 * x
        let (c0, c1) = poly1fit(
            &vec![1.0, 2.0, 3.0, 4.0, 5.0],
            &vec![4.0, 7.0, 10.0, 13.0, 16.0],
        );
        assert_eq!(c0, 1.0);
        assert_eq!(c1, 3.0);
        // test poly2fit -> y = 1 + 3 * x + 5 * x * x
        let (c0, c1, c2) = poly2fit(
            &vec![1.0, 2.0, 3.0, 4.0, 5.0],
            &vec![9.0, 27.0, 55.0, 93.0, 141.0],
        );
        assert_eq!(c0, 1.0);
        assert_eq!(c1, 3.0);
        assert_eq!(c2, 5.0);
    }
}
