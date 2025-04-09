"""HsTerm"""

import sympy as sp


z0, z1, z2, z3 = sp.symbols("z0,z1,z2,z3")
hs = (
    3 * z1 * z2 / (1 - z3) + z2**3 / z3 / (1 - z3) ** 2 + z2**3 / z3**2 * sp.log(1 - z3)
) / z0


hs_z0 = sp.diff(hs, z0)
print("$$ z0 =", sp.latex(hs_z0 * z0 * z0), "$$")


hs_z1 = sp.diff(hs, z1)
print("$$ z1 =", sp.latex(hs_z1 * z0), "$$")


hs_z2 = sp.diff(hs, z2)
print("$$ z2 =", sp.latex(hs_z2 * z0), "$$")


hs_z3 = sp.diff(hs, z3)
print("$$ z3 =", sp.latex(hs_z3 * z0), "$$")


hs_z = hs_z0 * z0 + hs_z1 * z1 + hs_z2 * z2 + hs_z3 * z3
print("$$ \\alpha^{hs}_\\rho =", sp.latex(hs_z.factor()), "$$")


hs_z0z0 = sp.diff(hs_z0, z0)
print("$$ z0z0 =", sp.latex(hs_z0z0 * z0 * z0 * z0), "$$")


hs_z0z1 = sp.diff(hs_z0, z1) + sp.diff(hs_z1, z0)
print("$$ z0z1 =", sp.latex(hs_z0z1 * z0 * z0), "$$")


hs_z0z2 = sp.diff(hs_z0, z2) + sp.diff(hs_z2, z0)
print("$$ z0z2 =", sp.latex(hs_z0z2 * z0 * z0), "$$")


hs_z0z3 = sp.diff(hs_z0, z3) + sp.diff(hs_z3, z0)
print("$$ z0z3 =", sp.latex(hs_z0z3 * z0 * z0), "$$")


hs_z1z1 = sp.diff(hs_z1, z1)
print("$$ z1z1 =", sp.latex(hs_z1z1 * z0), "$$")


hs_z1z2 = sp.diff(hs_z1, z2) + sp.diff(hs_z2, z1)
print("$$ z1z2 =", sp.latex(hs_z1z2 * z0), "$$")


hs_z1z3 = sp.diff(hs_z1, z3) + sp.diff(hs_z3, z1)
print("$$ z1z3 =", sp.latex(hs_z1z3 * z0), "$$")


hs_z2z2 = sp.diff(hs_z2, z2)
print("$$ z2z2 =", sp.latex(hs_z2z2 * z0), "$$")


hs_z2z3 = sp.diff(hs_z2, z3) + sp.diff(hs_z3, z2)
print("$$ z2z3 =", sp.latex(hs_z2z3 * z0), "$$")


hs_z3z3 = sp.diff(hs_z3, z3)
print("$$ z3z3 =", sp.latex(hs_z3z3 * z0), "$$")


hs_z0z = z0 * (hs_z0z0 * z0 + hs_z0z1 * z1 + hs_z0z2 * z2 + hs_z0z3 * z3)
hs_z1z = z1 * (hs_z1z1 * z1 + hs_z1z2 * z2 + hs_z1z3 * z3)
hs_z2z = z2 * (hs_z2z2 * z2 + hs_z2z3 * z3)
hs_z3z = z3 * (hs_z3z3 * z3)


hs_zz = hs_z0z + hs_z1z + hs_z2z + hs_z3z
print("$$ \\alpha^{hs}_{\\rho\\rho} =", sp.latex(hs_zz.factor()), "$$")


hs_z0z0z0 = sp.diff(hs_z0z0, z0)
print("$$ z0z0z0 =", sp.latex(hs_z0z0z0 * z0 * z0 * z0 * z0), "$$")


hs_z0z0z1 = sp.diff(hs_z0z0, z1) + sp.diff(hs_z0z1, z0)
print("$$ z0z0z1 =", sp.latex(hs_z0z0z1 * z0 * z0 * z0), "$$")


hs_z0z0z2 = sp.diff(hs_z0z0, z2) + sp.diff(hs_z0z2, z0)
print("$$ z0z0z2 =", sp.latex(hs_z0z0z2 * z0 * z0 * z0), "$$")


hs_z0z0z3 = sp.diff(hs_z0z0, z3) + sp.diff(hs_z0z3, z0)
print("$$ z0z0z3 =", sp.latex(hs_z0z0z3 * z0 * z0 * z0), "$$")


hs_z0z1z1 = sp.diff(hs_z0z1, z1) + sp.diff(hs_z1z1, z0)
print("$$ z0z1z1 =", sp.latex(hs_z0z1z1 * z0 * z0), "$$")


hs_z0z1z2 = sp.diff(hs_z0z1, z2) + sp.diff(hs_z0z2, z1) + sp.diff(hs_z1z2, z0)
print("$$ z0z1z2 =", sp.latex(hs_z0z1z2 * z0 * z0), "$$")


hs_z0z1z3 = sp.diff(hs_z0z1, z3) + sp.diff(hs_z0z3, z1) + sp.diff(hs_z1z3, z0)
print("$$ z0z1z3 =", sp.latex(hs_z0z1z3 * z0 * z0), "$$")


hs_z0z2z2 = sp.diff(hs_z0z2, z2) + sp.diff(hs_z2z2, z0)
print("$$ z0z2z2 =", sp.latex(hs_z0z2z2 * z0 * z0), "$$")


hs_z0z2z3 = sp.diff(hs_z0z2, z3) + sp.diff(hs_z0z3, z2) + sp.diff(hs_z2z3, z0)
print("$$ z0z2z3 =", sp.latex(hs_z0z2z3 * z0 * z0), "$$")


hs_z0z3z3 = sp.diff(hs_z0z3, z3) + sp.diff(hs_z3z3, z0)
print("$$ z0z3z3 =", sp.latex(hs_z0z3z3 * z0 * z0), "$$")


hs_z1z1z1 = sp.diff(hs_z1z1, z1)
print("$$ z1z1z1 =", sp.latex(hs_z1z1z1 * z0), "$$")


hs_z1z1z2 = sp.diff(hs_z1z1, z2) + sp.diff(hs_z1z2, z1)
print("$$ z1z1z2 =", sp.latex(hs_z1z1z2 * z0), "$$")


hs_z1z1z3 = sp.diff(hs_z1z1, z3) + sp.diff(hs_z1z3, z1)
print("$$ z1z1z3 =", sp.latex(hs_z1z1z3 * z0), "$$")


hs_z1z2z2 = sp.diff(hs_z1z2, z2) + sp.diff(hs_z2z2, z1)
print("$$ z1z2z2 =", sp.latex(hs_z1z2z2 * z0), "$$")


hs_z1z2z3 = sp.diff(hs_z1z2, z3) + sp.diff(hs_z1z3, z2) + sp.diff(hs_z2z3, z1)
print("$$ z1z2z3 =", sp.latex(hs_z1z2z3 * z0), "$$")


hs_z1z3z3 = sp.diff(hs_z1z3, z3) + sp.diff(hs_z3z3, z1)
print("$$ z1z3z3 =", sp.latex(hs_z1z3z3 * z0), "$$")


hs_z2z2z2 = sp.diff(hs_z2z2, z2)
print("$$ z2z2z2 =", sp.latex(hs_z2z2z2 * z0), "$$")


hs_z2z2z3 = sp.diff(hs_z2z2, z3) + sp.diff(hs_z2z3, z2)
print("$$ z2z2z3 =", sp.latex(hs_z2z2z3 * z0), "$$")


hs_z2z3z3 = sp.diff(hs_z2z3, z3) + sp.diff(hs_z3z3, z2)
print("$$ z2z3z3 =", sp.latex(hs_z2z3z3 * z0), "$$")


hs_z3z3z3 = sp.diff(hs_z3z3, z3)
print("$$ z3z3z3 =", sp.latex(hs_z3z3z3 * z0), "$$")


hs_z0z0z = z0 * z0 * (hs_z0z0z0 * z0 + hs_z0z0z1 * z1 + hs_z0z0z2 * z2 + hs_z0z0z3 * z3)
hs_z0z1z = z0 * z1 * (hs_z0z1z1 * z1 + hs_z0z1z2 * z2 + hs_z0z1z3 * z3)
hs_z0z2z = z0 * z2 * (hs_z0z2z2 * z2 + hs_z0z2z3 * z3)
hs_z0z3z = z0 * z3 * (hs_z0z3z3 * z3)
hs_z0z = hs_z0z0z + hs_z0z1z + hs_z0z2z + hs_z0z3z


hs_z1z1z = z1 * z1 * (hs_z1z1z1 * z1 + hs_z1z1z2 * z2 + hs_z1z1z3 * z3)
hs_z1z2z = z1 * z2 * (hs_z1z2z2 * z2 + hs_z1z2z3 * z3)
hs_z1z3z = z1 * z3 * (hs_z1z3z3 * z3)
hs_z1z = hs_z1z1z + hs_z1z2z + hs_z1z3z


hs_z2z2z = z2 * z2 * (hs_z2z2z2 * z2 + hs_z2z2z3 * z3)
hs_z2z3z = z2 * z3 * (hs_z2z3z3 * z3)
hs_z2z = hs_z2z2z + hs_z2z3z


hs_z3z3z = z3 * z3 * (hs_z3z3z3 * z3)
hs_z3z = hs_z3z3z


hs_zzz = hs_z0z + hs_z1z + hs_z2z + hs_z3z
print("$$ \\alpha^{hs}_{\\rho\\rho\\rho} =", sp.latex(hs_zzz.factor()), "$$")
