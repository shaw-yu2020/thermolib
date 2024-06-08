"""
plot ceos
"""

import math
import numpy as np
import matplotlib.pyplot as plt

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

fig, (ax1, ax2) = plt.subplots(2, 1)

t = []
p = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    vdw.t_flash(T)
    t.append(vdw.T())
    p.append(vdw.p())
    rho_v.append(vdw.rho_v())
    rho_l.append(vdw.rho_l())
t = np.array(t)
p_s = np.array(p)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t, p, color="red", label="Vdw")
ax2.plot(rho_v, p, color="red", label="Vdw")
ax2.plot(rho_l, p, color="red")

t = []
p = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    rk.t_flash(T)
    t.append(rk.T())
    p.append(rk.p())
    rho_v.append(rk.rho_v())
    rho_l.append(rk.rho_l())
t = np.array(t)
p = np.array(p)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t, p, color="yellow", label="Rk")
ax2.plot(rho_v, p, color="yellow", label="Rk")
ax2.plot(rho_l, p, color="yellow")

t = []
p = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    srk.t_flash(T)
    t.append(srk.T())
    p.append(srk.p())
    rho_v.append(srk.rho_v())
    rho_l.append(srk.rho_l())
t = np.array(t)
p = np.array(p)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t, p, color="blue", label="Srk")
ax2.plot(rho_v, p, color="blue", label="Srk")
ax2.plot(rho_l, p, color="blue")

t = []
p = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    pr.t_flash(T)
    t.append(pr.T())
    p.append(pr.p())
    rho_v.append(pr.rho_v())
    rho_l.append(pr.rho_l())
t = np.array(t)
p = np.array(p)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t, p, color="green", label="Pr")
ax2.plot(rho_v, p, color="green", label="Pr")
ax2.plot(rho_l, p, color="green")

ax1.set_xlabel("temperature K")
ax1.set_ylabel("pressure Pa")
ax1.legend()

ax2.set_xlabel("density kg/m3")
ax2.set_ylabel("pressure Pa")
ax2.legend()

plt.tight_layout()

plt.show()
