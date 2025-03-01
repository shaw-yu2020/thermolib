"""brent_zero"""

import sys
import numpy as np


EPS = np.sqrt(sys.float_info.epsilon)


def brent_zero(f, a, b):  # pylint: disable=R0912
    """brent_zero"""
    (fa, fb) = (f(a), f(b))
    (c, fc) = (a, fa)
    d = b - a
    e = b - a
    while True:
        if abs(fc) < abs(fb):
            (a, fa) = (b, fb)
            (b, fb) = (c, fc)
            (c, fc) = (a, fa)
        tol = 2.0 * EPS * abs(b) + 10 * EPS
        m = 0.5 * (c - b)
        if abs(m) <= tol or fb == 0.0:
            break
        if abs(e) < tol or abs(fa) <= abs(fb):
            e = m
            d = e
        else:
            s = fb / fa
            if a == c:
                p = 2.0 * m * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if 0.0 < p:
                q = -q
            else:
                p = -p
            s = e
            e = d
            if 2.0 * p < 3.0 * m * q - abs(tol * q) and p < abs(0.5 * s * q):
                d = p / q
            else:
                e = m
                d = e
        (a, fa) = (b, fb)
        if tol < abs(d):
            b = b + d
        elif 0.0 < m:
            b = b + tol
        else:
            b = b - tol
        fb = f(b)
        if (0.0 < fb and 0.0 < fc) or (fb <= 0.0 and fc <= 0.0):
            (c, fc) = (a, fa)
            e = b - a
            d = e
    return b


def main():
    """main"""
    print("f(x) = sin(x) - x / 2")
    print(f"x = {brent_zero(lambda x: np.sin(x) - x / 2, 1.0, 2.0):.10e}")
    print("f(x) = 2 * x - exp(-x)")
    print(f"x = {brent_zero(lambda x: 2 * x - np.exp(-x), 0.0, 1.0):.10e}")
    print("f(x) = x * exp(-x)")
    print(f"x = { brent_zero(lambda x: x * np.exp(-x), -1.0, 0.5):.10e}")
    print("f(x) = exp(x) - 1 / 100 / x / x")
    print(f"x = {brent_zero(lambda x: np.exp(x) - 1 / 100 / x / x, 0.0001, 20.0):.10e}")
    print("f(x) = (x + 3) * (x - 1) * (x - 1)")
    print(f"x = {brent_zero(lambda x: (x + 3) * (x - 1) * (x - 1), -5.0, 2.0):.10e}")


if __name__ == "__main__":
    main()
