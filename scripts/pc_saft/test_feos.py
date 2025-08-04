"""test_feos"""

import numpy as np
import si_units as si
from feos.pcsaft import PcSaftRecord  # pylint: disable=E0401,E0611
from feos.pcsaft import PcSaftParameters  # pylint: disable=E0401,E0611
from feos.eos import EquationOfState  # pylint: disable=E0401,E0611
from feos.eos import State  # pylint: disable=E0401,E0611
from feos.eos import PhaseEquilibrium  # pylint: disable=E0401,E0611


from thermolib import PcSaftMix2  # pylint:disable=no-name-in-module


def test_pc_saft_gly_pure():
    """test_pc_saft_gly_pure"""
    recoed = PcSaftRecord(
        1.5255, 3.23, 188.9, kappa_ab=0.035176, epsilon_k_ab=2899.5, na=1, nb=1
    )  # methanol
    parameters = PcSaftParameters.from_model_records([recoed])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 24676:
        print("Error in tp_flash :: rho() # methanol")


def test_pc_saft_pure():
    """test_pc_saft_pure"""
    recoed = PcSaftRecord(2.8611, 2.6826, 205.35)  # SO2
    parameters = PcSaftParameters.from_model_records([recoed])
    eos = EquationOfState.pcsaft(parameters)
    critical_point = State.critical_point(eos)
    if np.round(critical_point.temperature / si.KELVIN) != 438:
        print("Error in c_flash :: T() # SO2")
    if np.round(critical_point.pressure() / (1e3 * si.PASCAL)) != 9099:
        print("Error in c_flash :: p() # SO2")
    if np.round(critical_point.density / (si.MOL / si.METER**3)) != 8079:
        print("Error in c_flash :: rho() # SO2")
    pressure = PhaseEquilibrium.vapor_pressure(eos, 298.15 * si.KELVIN)
    if np.round(pressure[0] / si.PASCAL) != 396865:
        print("Error in t_flash :: p_s() # SO2")
    state_v = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=pressure[0],
        density_initialization="vapor",
    )
    if np.round(state_v.density / (si.MOL / si.METER**3)) != 169:
        print("Error in t_flash :: rho_v() # SO2")
    state_l = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=pressure[0],
        density_initialization="liquid",
    )
    if np.round(state_l.density / (si.MOL / si.METER**3)) != 21115:
        print("Error in t_flash :: rho_l() # SO2")
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 41:
        print("Error in tp_flash :: rho() # SO2")
    recoed = PcSaftRecord(
        1.0656, 3.0007, 366.51, kappa_ab=0.034868, epsilon_k_ab=2500.7, na=1, nb=1
    )  # H20
    parameters = PcSaftParameters.from_model_records([recoed])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 51179:
        print("Error in tp_flash :: rho() # H20")
    recoed = PcSaftRecord(1.5131, 3.1869, 163.33, q=4.4)  # CO2
    parameters = PcSaftParameters.from_model_records([recoed])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 41:
        print("Error in tp_flash :: rho() # CO2")
    recoed = PcSaftRecord(2.7447, 3.2742, 232.99, mu=2.88)  # Acetone
    parameters = PcSaftParameters.from_model_records([recoed])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 13337:
        print("Error in tp_flash :: rho() # Acetone")


def test_pc_saft_mix2():  # pylint: disable=too-many-branches
    """test_pc_saft_mix2"""
    methane = PcSaftRecord(1.0000, 3.7039, 150.03)  # methane (Gross and Sadowski 2001)
    ethane = PcSaftRecord(1.6069, 3.5206, 191.42)  # ethane (Gross and Sadowski 2001)
    parameters = PcSaftParameters.from_model_records([methane, ethane])
    eos = EquationOfState.pcsaft(parameters)
    fluids = PcSaftMix2(
        [0.5, 0.5], [1.0000, 1.6069], [3.7039, 3.5206], [150.03, 191.42], 0
    )
    temp = 199.92
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
    for xyp in xyp_exp:  # check tx_flash
        try:
            py = fluids.tx_flash(temp, xyp[0])
        except RuntimeError:
            ok_thermolib = False
        else:
            ok_thermolib = True
        try:
            bubble = PhaseEquilibrium.bubble_point(
                eos, temp * si.KELVIN, np.array([xyp[0], 1.0 - xyp[0]])
            )
        except RuntimeError:
            ok_feos = False
        else:
            ok_feos = True
        if ok_thermolib and ok_feos:
            if np.round(py[1] * 1e6) != np.round(bubble.vapor.molefracs[0] * 1e6):
                print(f"Error in tx_flash :: {py[1]:.6f} (thermolib::tx_flash)")
                print(
                    f"Error in tx_flash :: {bubble.vapor.molefracs[0]:.6f} (feos::bubble_point)"
                )
        elif ok_feos and ~ok_thermolib:
            print(f"x = {xyp[0]} :: Ok(bubble_point) but Err(tx_flash)")
    for xyp in xyp_exp:  # check ty_flash
        try:
            px = fluids.ty_flash(temp, xyp[1])
        except RuntimeError:
            ok_thermolib = False
        else:
            ok_thermolib = True
        try:
            dew = PhaseEquilibrium.dew_point(
                eos, temp * si.KELVIN, np.array([xyp[1], 1.0 - xyp[1]])
            )
        except RuntimeError:
            ok_feos = False
        else:
            ok_feos = True
        if ok_thermolib and ok_feos:
            if np.round(px[1] * 1e6) != np.round(dew.liquid.molefracs[0] * 1e6):
                print(f"Error in ty_flash :: {px[1]:.6f} (thermolib::ty_flash)")
                print(
                    f"Error in ty_flash :: {dew.liquid.molefracs[0]:.6f} (feos::dew_point)"
                )
        elif ok_feos and ~ok_thermolib:
            print(f"y = {xyp[0]} :: Ok(dew_point) but Err(ty_flash)")


def main():
    """main"""
    test_pc_saft_gly_pure()
    test_pc_saft_pure()
    test_pc_saft_mix2()


if __name__ == "__main__":
    main()
