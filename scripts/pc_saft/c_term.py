"""C_term"""

import sympy as sp


n = sp.symbols("n")
m = (8 * n - 2 * n**2) / (1 - n) ** 4
m1 = (20 * n - 27 * n**2 + 12 * n**3 - 2 * n**4) / (n**2 - 3 * n + 2) ** 2


m_n = sp.diff(m, n)
print("$m_{n}$")
print("$$", sp.latex(m_n.factor()), "$$")


m1_n = sp.diff(m1, n)
print("$m1_{n}$")
print("$$", sp.latex(m1_n.factor()), "$$")


m_nn = sp.diff(m_n, n)
print("$m_{nn}$")
print("$$", sp.latex(m_nn.factor()), "$$")


m1_nn = sp.diff(m1_n, n)
print("$m1_{nn}$")
print("$$", sp.latex(m1_nn.factor()), "$$")


m_nnn = sp.diff(m_nn, n)
print("$m_{nnn}$")
print("$$", sp.latex(m_nnn.factor()), "$$")


m1_nnn = sp.diff(m1_nn, n)
print("$m1_{nnn}$")
print("$$", sp.latex(m1_nnn.factor()), "$$")


m_nnnn = sp.diff(m_nnn, n)
print("$m_{nnnn}$")
print("$$", sp.latex(m_nnnn.factor()), "$$")


m1_nnnn = sp.diff(m1_nnn, n)
print("$m1_{nnnn}$")
print("$$", sp.latex(m1_nnnn.factor()), "$$")
