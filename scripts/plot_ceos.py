"""plot ceos"""

import math
import numpy as np
import matplotlib.pyplot as plt
from thermolib import hello, Vdw, Rk, Srk, Pr  # pylint: disable=no-name-in-module


print("test_hello:", hello())


TC = 430.64
PC = 7886600
OMEGA = 0.256
vdw = Vdw(TC, PC)
rk = Rk(TC, PC)
srk = Srk(TC, PC, OMEGA)
pr = Pr(TC, PC, OMEGA)


fig, (ax1, ax2) = plt.subplots(2, 1)


t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.6 * TC), math.ceil(TC)):
    vdw.t_flash(T)
    t_s.append(vdw.T_s())
    p_s.append(vdw.p_s())
    rho_v.append(vdw.rho_v())
    rho_l.append(vdw.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)
ax1.plot(t_s, p_s, color="r", label="Vdw")
ax2.plot(rho_v, p_s, color="r", label="Vdw")
ax2.plot(rho_l, p_s, color="r")


t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.6 * TC), math.ceil(TC)):
    rk.t_flash(T)
    t_s.append(rk.T_s())
    p_s.append(rk.p_s())
    rho_v.append(rk.rho_v())
    rho_l.append(rk.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)
ax1.plot(t_s, p_s, color="g", label="Rk")
ax2.plot(rho_v, p_s, color="g", label="Rk")
ax2.plot(rho_l, p_s, color="g")


t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.6 * TC), math.ceil(TC)):
    srk.t_flash(T)
    t_s.append(srk.T_s())
    p_s.append(srk.p_s())
    rho_v.append(srk.rho_v())
    rho_l.append(srk.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)
ax1.plot(t_s, p_s, color="b", label="Srk")
ax2.plot(rho_v, p_s, color="b", label="Srk")
ax2.plot(rho_l, p_s, color="b")


t_s = []
p_s = []
rho_v = []
rho_l = []
for T in range(math.floor(0.6 * TC), math.ceil(TC)):
    pr.t_flash(T)
    t_s.append(pr.T_s())
    p_s.append(pr.p_s())
    rho_v.append(pr.rho_v())
    rho_l.append(pr.rho_l())
t_s = np.array(t_s)
p_s = np.array(p_s)
rho_v = np.array(rho_v)
rho_l = np.array(rho_l)
ax1.plot(t_s, p_s, color="y", label="Pr")
ax2.plot(rho_v, p_s, color="y", label="Pr")
ax2.plot(rho_l, p_s, color="y")


ax1.set_xlabel("temperature K")
ax1.set_ylabel("pressure Pa")
ax1.legend()
ax2.set_xlabel("density mol/m3")
ax2.set_ylabel("pressure Pa")
ax2.legend()


plt.tight_layout()
plt.show()
