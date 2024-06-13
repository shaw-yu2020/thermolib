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

t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    vdw.t_flash(T)
    t_s.append(vdw.T_s())
    p_s.append(vdw.p_s())
    rho_v.append(vdw.rho_v())
    rho_l.append(vdw.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t_s, p_s, color="red", label="Vdw")
ax2.plot(rho_v, p_s, color="red", label="Vdw")
ax2.plot(rho_l, p_s, color="red")

t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    rk.t_flash(T)
    t_s.append(rk.T_s())
    p_s.append(rk.p_s())
    rho_v.append(rk.rho_v())
    rho_l.append(rk.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t_s, p_s, color="yellow", label="Rk")
ax2.plot(rho_v, p_s, color="yellow", label="Rk")
ax2.plot(rho_l, p_s, color="yellow")

t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    srk.t_flash(T)
    t_s.append(srk.T_s())
    p_s.append(srk.p_s())
    rho_v.append(srk.rho_v())
    rho_l.append(srk.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t_s, p_s, color="blue", label="Srk")
ax2.plot(rho_v, p_s, color="blue", label="Srk")
ax2.plot(rho_l, p_s, color="blue")

t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.7 * TC), math.ceil(TC)):
    pr.t_flash(T)
    t_s.append(pr.T_s())
    p_s.append(pr.p_s())
    rho_v.append(pr.rho_v())
    rho_l.append(pr.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)

ax1.plot(t_s, p_s, color="green", label="Pr")
ax2.plot(rho_v, p_s, color="green", label="Pr")
ax2.plot(rho_l, p_s, color="green")

ax1.set_xlabel("temperature K")
ax1.set_ylabel("pressure Pa")
ax1.legend()

ax2.set_xlabel("density kg/m3")
ax2.set_ylabel("pressure Pa")
ax2.legend()

plt.tight_layout()

plt.show()
