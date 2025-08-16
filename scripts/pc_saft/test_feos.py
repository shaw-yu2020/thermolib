"""test_feos"""

import numpy as np
import si_units as si  # pylint: disable=E0401


from feos.pcsaft import PcSaftRecord  # pylint: disable=E0401,E0611
from feos.pcsaft import PcSaftParameters  # pylint: disable=E0401,E0611
from feos.eos import EquationOfState  # pylint: disable=E0401,E0611
from feos.eos import State  # pylint: disable=E0401,E0611
from feos.eos import PhaseEquilibrium  # pylint: disable=E0401,E0611


SO2 = PcSaftRecord(2.8611, 2.6826, 205.35)  # sulfur_dioxide
H2O_2B = PcSaftRecord(
    1.0656, 3.0007, 366.51, kappa_ab=0.034868, epsilon_k_ab=2500.7, na=1, nb=1
)  # water
H2O_3B = PcSaftRecord(
    1.0656, 3.0007, 366.51, kappa_ab=0.034868, epsilon_k_ab=2500.7, na=1, nb=2
)  # water
H2O_4C = PcSaftRecord(
    1.0656, 3.0007, 366.51, kappa_ab=0.034868, epsilon_k_ab=2500.7, na=2, nb=2
)  # water
CO2_QQ = PcSaftRecord(1.5131, 3.1869, 163.33, q=4.4)  # carbon_dioxide
ACETONE_DD = PcSaftRecord(2.7447, 3.2742, 232.99, mu=2.88)  # acetone
CH3OH_2B = PcSaftRecord(
    1.5255, 3.2300, 188.90, kappa_ab=0.035176, epsilon_k_ab=2899.5, na=1, nb=1
)  # methanol
CH3OH_3B = PcSaftRecord(
    1.5255, 3.2300, 188.90, kappa_ab=0.035176, epsilon_k_ab=2899.5, na=1, nb=2
)  # methanol
CH4 = PcSaftRecord(1.0000, 3.7039, 150.03)  # methane (Gross and Sadowski 2001)
C2H6 = PcSaftRecord(1.6069, 3.5206, 191.42)  # ethane (Gross and Sadowski 2001)
CO2 = PcSaftRecord(2.0729, 2.7852, 169.21)  # carbon_dioxide
CL2_QQ = PcSaftRecord(1.4682, 3.4480, 269.67, q=3.0724)  # chlorine
C6H6_QQ = PcSaftRecord(2.2463, 3.7852, 296.24, q=5.5907)  # benzene
BUTANONE_DD = PcSaftRecord(2.9835, 3.4239, 244.99, mu=2.78)  # butanone
PROPANAL_DD = PcSaftRecord(2.6001, 3.2872, 235.21, mu=2.72)  # propanal


def test_pc_saft_pure():  # pylint: disable=too-many-statements,disable=too-many-branches
    """test_pc_saft_pure"""
    # SO2
    parameters = PcSaftParameters.from_model_records([SO2])
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
    # H2O_2B
    parameters = PcSaftParameters.from_model_records([H2O_2B])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 51179:
        print("Error in tp_flash :: rho() # H20_2B")
    # H2O_3B
    parameters = PcSaftParameters.from_model_records([H2O_3B])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 51371:
        print("Error in tp_flash :: rho() # H20_3B")
    # H2O_4C
    parameters = PcSaftParameters.from_model_records([H2O_4C])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 55670:
        print("Error in tp_flash :: rho() # H20_4C")
    # CO2_QQ
    parameters = PcSaftParameters.from_model_records([CO2_QQ])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 41:
        print("Error in tp_flash :: rho() # CO2")
    # ACETONE_DD
    parameters = PcSaftParameters.from_model_records([ACETONE_DD])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 13337:
        print("Error in tp_flash :: rho() # ACEONE_DD")
    # CH3OH_2B
    parameters = PcSaftParameters.from_model_records([CH3OH_2B])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([1.0]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 24676:
        print("Error in tp_flash :: rho() # CH3OH_2B")
    # CH3OH_3B
    parameters = PcSaftParameters.from_model_records([CH3OH_3B])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([1.0]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 24836:
        print("Error in tp_flash :: rho() # CH3OH_3B")


def test_pc_saft_mix2():  # pylint: disable=too-many-statements,disable=too-many-branches
    """test_pc_saft_mix2"""
    # CH4+C2H6
    parameters = PcSaftParameters.from_model_records([CH4, C2H6])
    eos = EquationOfState.pcsaft(parameters)
    temp = 199.92 * si.KELVIN
    # check tx_flash
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.0214, 0.9786]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.3258:
        print(f"Error in tx_flash :: x = {0.0214}, # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.0512, 0.9488]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.5351:
        print(f"Error in tx_flash :: x = {0.0512}, # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.1039, 0.8961]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.6986:
        print(f"Error in tx_flash :: x = {0.1039}, # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.1875, 0.8125]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.8054:
        print(f"Error in tx_flash :: x = {0.1875}, # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.3100, 0.6900]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.8712:
        print(f"Error in tx_flash :: x = {0.3100}, # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.4526, 0.5474]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.9079:
        print(f"Error in tx_flash :: x = {0.4526}, # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.6601, 0.3399]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.9371:
        print(f"Error in tx_flash :: x = {0.6601}, # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(eos, temp, np.array([0.7852, 0.2148]))
    if np.round(bubble.vapor.molefracs[0] * 1e4) / 1e4 != 0.9498:
        print(f"Error in tx_flash :: x = {0.7852}, # CH4+C2H6")
    # check ty_flash
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.3005, 0.6995]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.0190:
        print(f"Error in ty_flash :: y = {0.3005}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.5098, 0.4902]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.0462:
        print(f"Error in ty_flash :: y = {0.5098}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.6800, 0.3200]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.0951:
        print(f"Error in ty_flash :: y = {0.6800}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.7957, 0.2043]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.1763:
        print(f"Error in ty_flash :: y = {0.7957}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.8679, 0.1321]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.3008:
        print(f"Error in ty_flash :: y = {0.8679}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.9052, 0.0948]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.4384:
        print(f"Error in ty_flash :: y = {0.9052}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.9337, 0.0663]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.6296:
        print(f"Error in ty_flash :: y = {0.9337}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, temp, np.array([0.9461, 0.0539]))
    if np.round(dew.liquid.molefracs[0] * 1e4) / 1e4 != 0.7474:
        print(f"Error in ty_flash :: y = {0.9461}, # CH4+C2H6")
    # CH4+C2H6
    parameters = PcSaftParameters.from_model_records([CH4, C2H6])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 41:
        print("Error in tp_flash :: rho() # CH4+C2H6")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 198.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.92057:
        print(f"Error in tx_flash :: x = {0.5}, # CH4+C2H6")
    dew = PhaseEquilibrium.dew_point(eos, 198.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.04226:
        print(f"Error in ty_flash :: y = {0.5}, # CH4+C2H6")
    # CO2+CH3OH_2B
    parameters = PcSaftParameters.from_model_records([CO2, CH3OH_2B])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 45:
        print("Error in tp_flash :: rho() # CO2+CH3OH_2B")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.99649:
        print(f"Error in tx_flash :: x = {0.5}, # CO2+CH3OH_2B")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.00143:
        print(f"Error in ty_flash :: y = {0.5}, # CO2+CH3OH_2B")
    # CO2_QQ+CH3OH_2B
    parameters = PcSaftParameters.from_model_records([CO2_QQ, CH3OH_2B])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 45:
        print("Error in tp_flash :: rho() # CO2_QQ+CH3OH_2B")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.99624:
        print(f"Error in tx_flash :: x = {0.5}, # CO2_QQ+CH3OH_2B")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.00061:
        print(f"Error in ty_flash :: y = {0.5}, # CO2_QQ+CH3OH_2B")
    # ACETONE_DD+CH3OH_2B
    parameters = PcSaftParameters.from_model_records([ACETONE_DD, CH3OH_2B])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 17116:
        print("Error in tp_flash :: rho() # ACETONE_DD+CH3OH_2B")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.66555:
        print(f"Error in tx_flash :: x = {0.5}, # ACETONE_DD+CH3OH_2B")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.06152:
        print(f"Error in ty_flash :: y = {0.5}, # ACETONE_DD+CH3OH_2B")
    # Cl2_QQ+C6H6_QQ
    parameters = PcSaftParameters.from_model_records([CL2_QQ, C6H6_QQ])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 14187:  # true=14137
        print("Error in tp_flash :: rho() # Cl2_QQ+C6H6_QQ")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.97885:  # true=0.97915
        print(f"Error in tx_flash :: x = {0.5}, # Cl2_QQ+C6H6_QQ")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.02334:  # true=0.02098
        print(f"Error in ty_flash :: y = {0.5}, # Cl2_QQ+C6H6_QQ")
    # BUTANONE_DD+PROPANAL_DD
    parameters = PcSaftParameters.from_model_records([BUTANONE_DD, PROPANAL_DD])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 12165:
        print("Error in tp_flash :: rho() # BUTANONE_DD+PROPANAL_DD")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.22132:
        print(f"Error in tx_flash :: x = {0.5}, # BUTANONE_DD+PROPANAL_DD")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.78066:
        print(f"Error in ty_flash :: y = {0.5}, # BUTANONE_DD+PROPANAL_DD")
    # CO2_QQ+ACETONE_DD
    parameters = PcSaftParameters.from_model_records([CO2_QQ, ACETONE_DD])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.2, 0.8]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 14545:  # true=41?
        print("Error in tp_flash :: rho() # CO2_QQ+ACETONE_DD")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.98973:
        print(f"Error in tx_flash :: x = {0.5}, # CO2_QQ+ACETONE_DD")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.00688:
        print(f"Error in ty_flash :: y = {0.5}, # CO2_QQ+ACETONE_DD")


def main():
    """main"""
    test_pc_saft_pure()
    test_pc_saft_mix2()


if __name__ == "__main__":
    main()
