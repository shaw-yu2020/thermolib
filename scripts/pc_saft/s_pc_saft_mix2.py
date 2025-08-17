"""s_pc_saft_mix2"""

from thermolib import SPcSaftMix2  # pylint:disable=no-name-in-module


import numpy as np


import matplotlib.pyplot as plt


fluids = SPcSaftMix2(
    [0.5, 0.5], [1.0000, 1.6069], [3.7039, 3.5206], [150.03, 191.42], 0
)


TEMP = 199.92


xyp_exp = np.array(
    [
        [0.0214, 0.3005, 362e3],
        [0.0512, 0.5098, 442e3],
        [0.1039, 0.6800, 680e3],
        [0.1875, 0.7957, 1090e3],
        [0.3100, 0.8679, 1700e3],
        [0.4526, 0.9052, 2380e3],
        [0.6601, 0.9337, 3400e3],
        [0.7852, 0.9461, 4080e3],
        [0.8942, 0.9562, 4765e3],
        [0.9126, 0.9584, 4890e3],
        [0.9175, 0.9578, 4940e3],
        [0.9222, 0.9575, 4980e3],
        [0.9319, 0.9577, 5035e3],
    ]
)


xyp_tpz = []
for p in range(1, 61):
    try:
        xy = fluids.tpz_flash(TEMP, p * 1e5)
        xyp_tpz.append([xy[0], xy[1], p * 1e5])
    except RuntimeError:
        print(f"tpz_flash diverge in {p*1e5} Pa")
        continue
xyp_tpz = np.array(xyp_tpz)


xyp_tx = []
for x in range(0, 101):
    try:
        py = fluids.tx_flash(TEMP, x / 100)
        xyp_tx.append([x / 100, py[1], py[0]])
    except RuntimeError:
        print(f"tx_flash diverge in x = {x/100}")
        continue
xyp_tx = np.array(xyp_tx)


xyp_ty = []
for y in range(0, 101):
    try:
        px = fluids.ty_flash(TEMP, y / 100)
        xyp_ty.append([px[1], y / 100, px[0]])
    except RuntimeError:
        print(f"ty_flash diverge in y = {y/100}")
        continue
xyp_ty = np.array(xyp_ty)


plt.plot(xyp_tpz[:, 0], xyp_tpz[:, 2], c="r", label="tpz_flash")
plt.plot(xyp_tpz[:, 1], xyp_tpz[:, 2], c="r")


plt.scatter(xyp_tx[:, 0], xyp_tx[:, 2], c="g", label="tx_flash")
plt.scatter(xyp_tx[:, 1], xyp_tx[:, 2], c="g")


plt.scatter(xyp_ty[:, 0], xyp_ty[:, 2], c="b", label="ty_flash")
plt.scatter(xyp_ty[:, 1], xyp_ty[:, 2], c="b")


plt.scatter(xyp_exp[:, 0], xyp_exp[:, 2], c="w", edgecolors="k", marker="o")
plt.scatter(xyp_exp[:, 1], xyp_exp[:, 2], c="w", edgecolors="k", marker="o")


plt.legend()
plt.show()
