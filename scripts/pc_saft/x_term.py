"""X_term"""

import sympy as sp


X = sp.symbols("X")
b = sp.symbols("b")


# AssocType::2B => XA and XB
t = (1 - X) / X**2


tX = sp.diff(t, X)
print("$t_{X}$")
print("$$", sp.latex(tX.factor()), "$$")


Xt = 1 / tX
print("$X_{t}$")
print("$$", sp.latex(Xt.factor()), "$$")


Xtt = sp.diff(Xt, X) * Xt
print("$X_{tt}$")
print("$$", sp.latex(Xtt.factor()), "$$")


Xttt = sp.diff(Xtt, X) * Xt
print("$X_{ttt}$")
print("$$", sp.latex(Xttt.factor()), "$$")


Xtttt = sp.diff(Xttt, X) * Xt
print("$X_{tttt}$")
print("$$", sp.latex(Xtttt.factor()), "$$")


# AssocType::3B => XA and XB
t = (1 - X) / (2 * X**2 - X)


tX = sp.diff(t, X)
print("$t_{X}$")
print("$$", sp.latex(tX.factor()), "$$")


Xt = 1 / tX
print("$X_{t}$")
print("$$", sp.latex(Xt.factor()), "$$")


Xtt = sp.diff(Xt, X) * Xt
print("$X_{tt}$")
print("$$", sp.latex(Xtt.factor()), "$$")


Xttt = sp.diff(Xtt, X) * Xt
print("$X_{ttt}$")
print("$$", sp.latex(Xttt.factor()), "$$")


Xtttt = sp.diff(Xttt, X) * Xt
print("$X_{tttt}$")
print("$$", sp.latex(Xtttt.factor()), "$$")


# AssocType::3B => XC
t = (1 - X) / (X**2 + X)


tX = sp.diff(t, X)
print("$t_{X}$")
print("$$", sp.latex(tX.factor()), "$$")


Xt = 1 / tX
print("$X_{t}$")
print("$$", sp.latex(Xt.factor()), "$$")


Xtt = sp.diff(Xt, X) * Xt
print("$X_{tt}$")
print("$$", sp.latex(Xtt.factor()), "$$")


Xttt = sp.diff(Xtt, X) * Xt
print("$X_{ttt}$")
print("$$", sp.latex(Xttt.factor()), "$$")


Xtttt = sp.diff(Xttt, X) * Xt
print("$X_{tttt}$")
print("$$", sp.latex(Xtttt.factor()), "$$")


# AssocType::4C => XA and XB and XC and XD
t = (1 - X) / X**2 / 2


tX = sp.diff(t, X)
print("$t_{X}$")
print("$$", sp.latex(tX.factor()), "$$")


Xt = 1 / tX
print("$X_{t}$")
print("$$", sp.latex(Xt.factor()), "$$")


Xtt = sp.diff(Xt, X) * Xt
print("$X_{tt}$")
print("$$", sp.latex(Xtt.factor()), "$$")


Xttt = sp.diff(Xtt, X) * Xt
print("$X_{ttt}$")
print("$$", sp.latex(Xttt.factor()), "$$")


Xtttt = sp.diff(Xttt, X) * Xt
print("$X_{tttt}$")
print("$$", sp.latex(Xtttt.factor()), "$$")


# AssocType::3B+ => XA
t = (1 - X) / ((1 + b) * X**2 - b * X)


tX = sp.diff(t, X)
print("$t_{X}$")
print("$$", sp.latex(tX.factor()), "$$")


Xt = 1 / tX
print("$X_{t}$")
print("$$", sp.latex(Xt.factor()), "$$")


Xtt = sp.diff(Xt, X) * Xt
print("$X_{tt}$")
print("$$", sp.latex(Xtt.factor()), "$$")
XttC = 2 * X**3 * (X * b + X - b) ** 3 / (X**2 * b + X**2 - 2 * X * b - 2 * X + b) ** 3
XttAll = (Xtt / XttC).factor(X)
print("$$", sp.latex(XttAll.simplify()), "$$")
print("$$X^3=", sp.latex(XttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XttAll.coeff(X**1).factor()), "$$")


Xttt = sp.diff(Xtt, X) * Xt
print("$X_{ttt}$")
print("$$", sp.latex(Xttt.factor()), "$$")
XtttC = 6 * X**4 * (X * b + X - b) ** 4 / (X**2 * b + X**2 - 2 * X * b - 2 * X + b) ** 5
XtttAll = (Xttt / XtttC).factor(X)
print("$$", sp.latex(XtttAll.simplify()), "$$")
print("$$X^6=", sp.latex(XtttAll.coeff(X**6).factor()), "$$")
print("$$X^5=", sp.latex(XtttAll.coeff(X**5).factor()), "$$")
print("$$X^4=", sp.latex(XtttAll.coeff(X**4).factor()), "$$")
print("$$X^3=", sp.latex(XtttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XtttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XtttAll.coeff(X**1).factor()), "$$")


Xtttt = sp.diff(Xttt, X) * Xt
print("$X_{tttt}$")
print("$$", sp.latex(Xtttt.factor()), "$$")
XttttC = (
    24 * X**5 * (X * b + X - b) ** 5 / (X**2 * b + X**2 - 2 * X * b - 2 * X + b) ** 7
)
XttttAll = (Xtttt / XttttC).factor(X)
print("$$", sp.latex(XttttAll.simplify()), "$$")
print("$$X^9=", sp.latex(XttttAll.coeff(X**9).factor()), "$$")
print("$$X^8=", sp.latex(XttttAll.coeff(X**8).factor()), "$$")
print("$$X^7=", sp.latex(XttttAll.coeff(X**7).factor()), "$$")
print("$$X^6=", sp.latex(XttttAll.coeff(X**6).factor()), "$$")
print("$$X^5=", sp.latex(XttttAll.coeff(X**5).factor()), "$$")
print("$$X^4=", sp.latex(XttttAll.coeff(X**4).factor()), "$$")
print("$$X^3=", sp.latex(XttttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XttttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XttttAll.coeff(X**1).factor()), "$$")


# AssocType::3B+ => XB
t = (b * (1 - X)) / ((1 + b) * X**2 + (b**2 - 2) * X + (1 - b))


tX = sp.diff(t, X)
print("$t_{X}$")
print("$$", sp.latex(tX.factor()), "$$")


Xt = 1 / tX
print("$X_{t}$")
print("$$", sp.latex(Xt.factor()), "$$")


Xtt = sp.diff(Xt, X) * Xt
print("$X_{tt}$")
print("$$", sp.latex(Xtt.factor()), "$$")
XttC = (
    -2
    * (X + b - 1) ** 3
    * (X * b + X - 1) ** 3
    / b**2
    / (-(X**2) * b - X**2 + 2 * X * b + 2 * X + b**2 - b - 1) ** 3
)
XttAll = (Xtt / XttC).factor(X)
print("$$", sp.latex(XttAll.simplify()), "$$")
print("$$X^3=", sp.latex(XttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XttAll.coeff(X**1).factor()), "$$")


Xttt = sp.diff(Xtt, X) * Xt
print("$X_{ttt}$")
print("$$", sp.latex(Xttt.factor()), "$$")
XtttC = (
    -6
    * (X + b - 1) ** 4
    * (X * b + X - 1) ** 4
    / b**3
    / (-(X**2) * b - X**2 + 2 * X * b + 2 * X + b**2 - b - 1) ** 5
)
XtttAll = (Xttt / XtttC).factor(X)
print("$$", sp.latex(XtttAll.simplify()), "$$")
print("$$X^6=", sp.latex(XtttAll.coeff(X**6).factor()), "$$")
print("$$X^5=", sp.latex(XtttAll.coeff(X**5).factor()), "$$")
print("$$X^4=", sp.latex(XtttAll.coeff(X**4).factor()), "$$")
print("$$X^3=", sp.latex(XtttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XtttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XtttAll.coeff(X**1).factor()), "$$")


Xtttt = sp.diff(Xttt, X) * Xt
print("$X_{tttt}$")
print("$$", sp.latex(Xtttt.factor()), "$$")
XttttC = (
    -24
    * (X + b - 1) ** 5
    * (X * b + X - 1) ** 5
    / b**4
    / (-(X**2) * b - X**2 + 2 * X * b + 2 * X + b**2 - b - 1) ** 7
)
XttttAll = (Xtttt / XttttC).factor(X)
print("$$", sp.latex(XttttAll.simplify()), "$$")
print("$$X^9=", sp.latex(XttttAll.coeff(X**9).factor()), "$$")
print("$$X^8=", sp.latex(XttttAll.coeff(X**8).factor()), "$$")
print("$$X^7=", sp.latex(XttttAll.coeff(X**7).factor()), "$$")
print("$$X^6=", sp.latex(XttttAll.coeff(X**6).factor()), "$$")
print("$$X^5=", sp.latex(XttttAll.coeff(X**5).factor()), "$$")
print("$$X^4=", sp.latex(XttttAll.coeff(X**4).factor()), "$$")
print("$$X^3=", sp.latex(XttttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XttttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XttttAll.coeff(X**1).factor()), "$$")


# AssocType::3B+ => XC
t = (1 - X) / (X**2 + b * X)


tX = sp.diff(t, X)
print("$t_{X}$")
print("$$", sp.latex(tX.factor()), "$$")


Xt = 1 / tX
print("$X_{t}$")
print("$$", sp.latex(Xt.factor()), "$$")


Xtt = sp.diff(Xt, X) * Xt
print("$X_{tt}$")
print("$$", sp.latex(Xtt.factor()), "$$")
XttC = -2 * X**3 * (X + b) ** 3 / (-(X**2) + 2 * X + b) ** 3
XttAll = (Xtt / XttC).factor(X)
print("$$", sp.latex(XttAll.simplify()), "$$")
print("$$X^3=", sp.latex(XttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XttAll.coeff(X**1).factor()), "$$")


Xttt = sp.diff(Xtt, X) * Xt
print("$X_{ttt}$")
print("$$", sp.latex(Xttt.factor()), "$$")
XtttC = -6 * X**4 * (X + b) ** 4 / (-(X**2) + 2 * X + b) ** 5
XtttAll = (Xttt / XtttC).factor(X)
print("$$", sp.latex(XtttAll.simplify()), "$$")
print("$$X^6=", sp.latex(XtttAll.coeff(X**6).factor()), "$$")
print("$$X^5=", sp.latex(XtttAll.coeff(X**5).factor()), "$$")
print("$$X^4=", sp.latex(XtttAll.coeff(X**4).factor()), "$$")
print("$$X^3=", sp.latex(XtttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XtttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XtttAll.coeff(X**1).factor()), "$$")


Xtttt = sp.diff(Xttt, X) * Xt
print("$X_{tttt}$")
print("$$", sp.latex(Xtttt.factor()), "$$")
XttttC = -24 * X**5 * (X + b) ** 5 / (-(X**2) + 2 * X + b) ** 7
XttttAll = (Xtttt / XttttC).factor(X)
print("$$", sp.latex(XttttAll.simplify()), "$$")
print("$$X^9=", sp.latex(XttttAll.coeff(X**9).factor()), "$$")
print("$$X^8=", sp.latex(XttttAll.coeff(X**8).factor()), "$$")
print("$$X^7=", sp.latex(XttttAll.coeff(X**7).factor()), "$$")
print("$$X^6=", sp.latex(XttttAll.coeff(X**6).factor()), "$$")
print("$$X^5=", sp.latex(XttttAll.coeff(X**5).factor()), "$$")
print("$$X^4=", sp.latex(XttttAll.coeff(X**4).factor()), "$$")
print("$$X^3=", sp.latex(XttttAll.coeff(X**3).factor()), "$$")
print("$$X^2=", sp.latex(XttttAll.coeff(X**2).factor()), "$$")
print("$$X^1=", sp.latex(XttttAll.coeff(X**1).factor()), "$$")
