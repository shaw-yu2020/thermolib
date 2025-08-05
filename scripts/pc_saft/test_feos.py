"""test_feos"""

import numpy as np
import si_units as si
from feos.pcsaft import PcSaftRecord  # pylint: disable=E0401,E0611
from feos.pcsaft import PcSaftParameters  # pylint: disable=E0401,E0611
from feos.eos import EquationOfState  # pylint: disable=E0401,E0611
from feos.eos import State  # pylint: disable=E0401,E0611
from feos.eos import PhaseEquilibrium  # pylint: disable=E0401,E0611


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


def test_pc_saft_mix2():  # pylint: disable=too-many-statements
    """test_pc_saft_mix2"""
    methane = PcSaftRecord(1.0000, 3.7039, 150.03)  # methane (Gross and Sadowski 2001)
    ethane = PcSaftRecord(1.6069, 3.5206, 191.42)  # ethane (Gross and Sadowski 2001)
    parameters = PcSaftParameters.from_model_records([methane, ethane])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 41:
        print("Error in tp_flash :: rho() # C1-C2")
    temp = 199.92 * si.KELVIN
    # check tx_flash
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.0214, 0.9786]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.3258:
        print(f"Error in tx_flash :: x = {0.0214}, # C1-C2")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.0512, 0.9488]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.5351:
        print(f"Error in tx_flash :: x = {0.0512}, # C1-C2")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.1039, 0.8961]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.6986:
        print(f"Error in tx_flash :: x = {0.1039}, # C1-C2")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.1875, 0.8125]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.8054:
        print(f"Error in tx_flash :: x = {0.1875}, # C1-C2")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.3100, 0.6900]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.8712:
        print(f"Error in tx_flash :: x = {0.3100}, # C1-C2")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.4526, 0.5474]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.9079:
        print(f"Error in tx_flash :: x = {0.4526}, # C1-C2")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.6601, 0.3399]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.9371:
        print(f"Error in tx_flash :: x = {0.6601}, # C1-C2")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.7852, 0.2148]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.9498:
        print(f"Error in tx_flash :: x = {0.7852}, # C1-C2")
    # check ty_flash
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.3005, 0.6995]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.0190:
        print(f"Error in ty_flash :: y = {0.3005}, # C1-C2")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.5098, 0.4902]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.0462:
        print(f"Error in ty_flash :: y = {0.5098}, # C1-C2")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.6800, 0.3200]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.0951:
        print(f"Error in ty_flash :: y = {0.6800}, # C1-C2")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.7957, 0.2043]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.1763:
        print(f"Error in ty_flash :: y = {0.7957}, # C1-C2")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.8679, 0.1321]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.3008:
        print(f"Error in ty_flash :: y = {0.8679}, # C1-C2")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.9052, 0.0948]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.4384:
        print(f"Error in ty_flash :: y = {0.9052}, # C1-C2")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.9337, 0.0663]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.6296:
        print(f"Error in ty_flash :: y = {0.9337}, # C1-C2")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.9461, 0.0539]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.7474:
        print(f"Error in ty_flash :: y = {0.9461}, # C1-C2")
    # test assoc_term
    carbon_dioxide = PcSaftRecord(
        2.0729, 2.7852, 169.21
    )  # carbon_dioxide (Gross and Sadowski 2001)
    methanol = PcSaftRecord(
        1.5255, 3.2300, 188.90, kappa_ab=0.035176, epsilon_k_ab=2899.5, na=1, nb=1
    )  # methanol (Gross and Sadowski 2002)
    parameters = PcSaftParameters.from_model_records([carbon_dioxide, methanol])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 45:
        print("Error in tp_flash :: rho() # carbon_dioxide + methanol")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.9965:
        print(f"Error in tx_flash :: x = {0.5}, # carbon_dioxide + methanol")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.0014:
        print(f"Error in ty_flash :: y = {0.5}, # carbon_dioxide + methanol")


def main():
    """main"""
    test_pc_saft_gly_pure()
    test_pc_saft_pure()
    test_pc_saft_mix2()


if __name__ == "__main__":
    main()
