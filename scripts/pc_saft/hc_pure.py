"""
HcPure
"""

import sympy as sp


eta = sp.symbols("eta")


hs = eta * (4 - 3 * eta) / (1 - eta) ** 2
print("$$ hs0 =", sp.latex(hs.factor()), "$$")


hs1 = sp.diff(hs, eta)
print("$$ hs1 =", sp.latex(hs1.factor()), "$$")


hs2 = sp.diff(hs1, eta)
print("$$ hs2 =", sp.latex(hs2.factor()), "$$")


hs3 = sp.diff(hs2, eta)
print("$$ hs3 =", sp.latex(hs3.factor()), "$$")


hs4 = sp.diff(hs3, eta)
print("$$ hs4 =", sp.latex(hs4.factor()), "$$")


gii = (2 - eta) / 2 / (1 - eta) ** 3
print("$$ gii0 =", sp.latex(gii.factor()), "$$")


gii1 = sp.diff(gii, eta)
print("$$ gii1 =", sp.latex(gii1.factor()), "$$")


gii2 = sp.diff(gii1, eta)
print("$$ gii2 =", sp.latex(gii2.factor()), "$$")


gii3 = sp.diff(gii2, eta)
print("$$ gii3 =", sp.latex(gii3.factor()), "$$")


gii4 = sp.diff(gii3, eta)
print("$$ gii4 =", sp.latex(gii4.factor()), "$$")


lng = sp.ln(gii)
print("$$ lng0 =", sp.latex(lng.simplify()), "$$")


lng1 = sp.diff(lng, eta)
print("$$ lng1 =", sp.latex(lng1.factor()), "$$")


lng2 = sp.diff(lng1, eta)
print("$$ lng2 =", sp.latex(lng2.factor()), "$$")


lng3 = sp.diff(lng2, eta)
print("$$ lng3 =", sp.latex(lng3.factor()), "$$")


lng4 = sp.diff(lng3, eta)
print("$$ lng4 =", sp.latex(lng4.factor()), "$$")
