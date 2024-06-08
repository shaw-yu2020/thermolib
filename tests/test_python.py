"""
test thermolib
"""

import math
import matplotlib.pyplot as plt
import numpy as np

# pylint: disable=E
from thermolib import hello, Vdw, Rk, Srk, Pr

print("test_hello:", hello())

TC = 430.64
PC = 7886600
OMEGA = 0.256
M = 0.064064

vdw = Vdw(TC, PC, M)
vdw.set_mass_unit()

rk = Rk(TC, PC, M)
rk.set_mass_unit()

srk = Srk(TC, PC, OMEGA, M)
srk.set_mass_unit()

pr = Pr(TC, PC, OMEGA, M)
pr.set_mass_unit()

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

p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    srk.t_flash(T)
    p_s.append(srk.p())
    rho_v.append(srk.rho_v())
    rho_l.append(srk.rho_l())
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

plt.plot(rho_v, p_s, color="blue")
plt.plot(rho_l, p_s, color="blue")

p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    pr.t_flash(T)
    p_s.append(pr.p())
    rho_v.append(pr.rho_v())
    rho_l.append(pr.rho_l())
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

plt.plot(rho_v, p_s, color="yellow")
plt.plot(rho_l, p_s, color="yellow")

plt.show()
