"""X_term"""

import sympy as sp


X = sp.symbols("X")


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
