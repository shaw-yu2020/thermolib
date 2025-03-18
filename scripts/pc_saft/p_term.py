"""PolarTerm"""

import sympy as sp


a2, a3 = sp.symbols("a_2 a_3")
a = a2**2 / (a2 - a3)


a2_d, a3_d = sp.symbols("a_{2d} a_{3d}")
a_d = sp.diff(a, a2) * a2_d + sp.diff(a, a3) * a3_d
print("$\\frac{a}{d}$")
print("$$", sp.latex(a_d.factor()), "$$")


a2_dd, a3_dd = sp.symbols("a_{2dd} a_{3dd}")
a_dd = (
    sp.diff(a_d, a2) * a2_d
    + sp.diff(a_d, a3) * a3_d
    + sp.diff(a_d, a2_d) * a2_dd
    + sp.diff(a_d, a3_d) * a3_dd
)
print("$\\frac{a}{dd}$")
print("$$", sp.latex(a_dd.factor()), "$$")


a2_ddd, a3_ddd = sp.symbols("a_{2ddd} a_{3ddd}")
a_ddd = (
    sp.diff(a_dd, a2) * a2_d
    + sp.diff(a_dd, a3) * a3_d
    + sp.diff(a_dd, a2_d) * a2_dd
    + sp.diff(a_dd, a3_d) * a3_dd
    + sp.diff(a_dd, a2_dd) * a2_ddd
    + sp.diff(a_dd, a3_dd) * a3_ddd
)
print("$\\frac{a}{ddd}$")
print("$$", sp.latex(a_ddd.factor()), "$$")


a2_dddd, a3_dddd = sp.symbols("a_{2dddd} a_{3dddd}")
a_dddd = (
    sp.diff(a_ddd, a2) * a2_d
    + sp.diff(a_ddd, a3) * a3_d
    + sp.diff(a_ddd, a2_d) * a2_dd
    + sp.diff(a_ddd, a3_d) * a3_dd
    + sp.diff(a_ddd, a2_dd) * a2_ddd
    + sp.diff(a_ddd, a3_dd) * a3_ddd
    + sp.diff(a_ddd, a2_ddd) * a2_dddd
    + sp.diff(a_ddd, a3_ddd) * a3_dddd
)
print("$\\frac{a}{dddd}$")
print("$$", sp.latex(a_dddd.factor()), "$$")


a2_t, a3_t = sp.symbols("a_{2t} a_{3t}")
a_t = sp.diff(a, a2) * a2_t + sp.diff(a, a3) * a3_t
print("$\\frac{a}{t}$")
print("$$", sp.latex(a_t.factor()), "$$")


a2_td, a3_td = sp.symbols("a_{2td} a_{3td}")
a_td = (
    sp.diff(a_t, a2) * a2_d
    + sp.diff(a_t, a3) * a3_d
    + sp.diff(a_t, a2_t) * a2_td
    + sp.diff(a_t, a3_t) * a3_td
)
print("$\\frac{a}{td}$")
print("$$", sp.latex(a_td.factor()), "$$")
a_dt = (
    sp.diff(a_d, a2) * a2_t
    + sp.diff(a_d, a3) * a3_t
    + sp.diff(a_d, a2_d) * a2_td
    + sp.diff(a_d, a3_d) * a3_td
)
print("$\\frac{a}{dt}$")
print("$$", sp.latex(a_dt.factor()), "$$")


a2_tdd, a3_tdd = sp.symbols("a_{2tdd} a_{3tdd}")
a_tdd = (
    sp.diff(a_td, a2) * a2_d
    + sp.diff(a_td, a3) * a3_d
    + sp.diff(a_td, a2_t) * a2_td
    + sp.diff(a_td, a3_t) * a3_td
    + sp.diff(a_td, a2_d) * a2_dd
    + sp.diff(a_td, a3_d) * a3_dd
    + sp.diff(a_td, a2_td) * a2_tdd
    + sp.diff(a_td, a3_td) * a3_tdd
)
print("$\\frac{a}{tdd}$")
print("$$", sp.latex(a_tdd.factor()), "$$")
a_ddt = (
    sp.diff(a_dd, a2) * a2_t
    + sp.diff(a_dd, a3) * a3_t
    + sp.diff(a_dd, a2_d) * a2_td
    + sp.diff(a_dd, a3_d) * a3_td
    + sp.diff(a_dd, a2_dd) * a2_tdd
    + sp.diff(a_dd, a3_dd) * a3_tdd
)
print("$\\frac{a}{ddt}$")
print("$$", sp.latex(a_ddt.factor()), "$$")


a2_tddd, a3_tddd = sp.symbols("a_{2tddd} a_{3tddd}")
a_tddd = (
    sp.diff(a_tdd, a2) * a2_d
    + sp.diff(a_tdd, a3) * a3_d
    + sp.diff(a_tdd, a2_t) * a2_td
    + sp.diff(a_tdd, a3_t) * a3_td
    + sp.diff(a_tdd, a2_d) * a2_dd
    + sp.diff(a_tdd, a3_d) * a3_dd
    + sp.diff(a_tdd, a2_td) * a2_tdd
    + sp.diff(a_tdd, a3_td) * a3_tdd
    + sp.diff(a_tdd, a2_dd) * a2_ddd
    + sp.diff(a_tdd, a3_dd) * a3_ddd
    + sp.diff(a_tdd, a2_tdd) * a2_tddd
    + sp.diff(a_tdd, a3_tdd) * a3_tddd
)
print("$\\frac{a}{tddd}$")
print("$$", sp.latex(a_tddd.factor()), "$$")
a_dddt = (
    sp.diff(a_ddd, a2) * a2_t
    + sp.diff(a_ddd, a3) * a3_t
    + sp.diff(a_ddd, a2_d) * a2_td
    + sp.diff(a_ddd, a3_d) * a3_td
    + sp.diff(a_ddd, a2_dd) * a2_tdd
    + sp.diff(a_ddd, a3_dd) * a3_tdd
    + sp.diff(a_ddd, a2_ddd) * a2_tddd
    + sp.diff(a_ddd, a3_ddd) * a3_tddd
)
print("$\\frac{a}{dddt}$")
print("$$", sp.latex(a_dddt.factor()), "$$")


a2_tt, a3_tt = sp.symbols("a_{2tt} a_{3tt}")
a_tt = (
    sp.diff(a_t, a2) * a2_t
    + sp.diff(a_t, a3) * a3_t
    + sp.diff(a_t, a2_t) * a2_tt
    + sp.diff(a_t, a3_t) * a3_tt
)
print("$\\frac{a}{tt}$")
print("$$", sp.latex(a_tt.factor()), "$$")
