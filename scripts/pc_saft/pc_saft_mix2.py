"""pc_saft_mix2"""

from feos.pcsaft import PcSaftRecord  # pylint: disable=E0401,E0611
from feos.pcsaft import PcSaftParameters  # pylint: disable=E0401,E0611
from feos.eos import EquationOfState  # pylint: disable=E0401,E0611
from feos.eos import PhaseEquilibrium  # pylint: disable=E0401,E0611


from thermolib import PcSaftMix2  # pylint:disable=no-name-in-module


import numpy as np
import si_units as si
import matplotlib.pyplot as plt


plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


# feos
CO2_QQ = PcSaftRecord(1.5131, 3.1869, 163.33, q=4.4)  # carbon_dioxide
ACETONE_DD = PcSaftRecord(2.7447, 3.2742, 232.99, mu=2.88)  # acetone
parameters = PcSaftParameters.from_model_records([CO2_QQ, ACETONE_DD])
eos = EquationOfState.pcsaft(parameters)


# thermolib
fluids = PcSaftMix2([0.5, 0.5], [1.5131, 2.7447], [3.1869, 3.2742], [163.33, 232.99], 0)
fluids.set_polar_term([4.4, 2.88], [4, 2])


# REF: https://doi.org/10.1016/j.supflu.2015.02.003
temp = np.array([313.15, 333.15, 353.15])
frac = np.array(range(47, 87, 1)) / 100
tpx = np.array(
    [
        [353.15, 10530, 0.794],
        [353.15, 9730, 0.729],
        [353.15, 8600, 0.648],
        [353.15, 7160, 0.552],
        [353.15, 6200, 0.485],
        [333.15, 8800, 0.823],
        [333.15, 7980, 0.756],
        [333.15, 7040, 0.675],
        [333.15, 5900, 0.582],
        [333.15, 5130, 0.516],
        [313.15, 6670, 0.843],
        [313.15, 6040, 0.778],
        [313.15, 5360, 0.699],
        [313.15, 4550, 0.612],
        [313.15, 4060, 0.551],
    ]
)


# feos::bubble_point
tpx_tx = []
for t in temp:
    for x in frac:
        try:
            bubble = PhaseEquilibrium.bubble_point(
                eos, t * si.KELVIN, np.array([x, 1 - x])
            )
            tpx_tx.append([t, bubble.vapor.pressure() / si.PASCAL, x])
        except RuntimeError:
            print("feos diverge")
            continue
tpx_tx = np.array(tpx_tx)
plt.scatter(tpx_tx[:, 2], tpx_tx[:, 1] / 1e3, c="r", marker="o", label="feos")


# thermolib::tx_flash
tpx_tx = []
for t in temp:
    for x in frac:
        try:
            py = fluids.tx_flash(t, x)
            tpx_tx.append([t, py[0], x])
        except RuntimeError:
            print("thermolib diverge")
            continue
tpx_tx = np.array(tpx_tx)
plt.scatter(tpx_tx[:, 2], tpx_tx[:, 1] / 1e3, c="g", marker="+", label="thermolib")


plt.scatter(tpx[:, 2], tpx[:, 1], c="w", edgecolors="k", marker="o")
plt.xlim((0.45, 0.88))
plt.xlabel("x(CO2)")
plt.ylabel("P(kPa)")
plt.legend()
plt.show()
