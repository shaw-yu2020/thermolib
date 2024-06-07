"""
test thermolib
"""

import math
import matplotlib.pyplot as plt
import numpy as np

# pylint: disable=E
from thermolib import hello, Vdw, Rk

print("test_hello:", hello())

TC = 430.64
PC = 7886600
M = 0.064064

vdw = Vdw(TC, PC, M)
vdw.set_mass_unit()

rk = Rk(TC, PC, M)
rk.set_mass_unit()

p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    vdw.t_flash(T)
    p_s.append(vdw.p())
    rho_v.append(vdw.rho_v())
    rho_l.append(vdw.rho_l())
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

plt.plot(rho_v, p_s, color="red")
plt.plot(rho_l, p_s, color="red")

p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    rk.t_flash(T)
    p_s.append(rk.p())
    rho_v.append(rk.rho_v())
    rho_l.append(rk.rho_l())
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

plt.plot(rho_v, p_s, color="green")
plt.plot(rho_l, p_s, color="green")

plt.show()
