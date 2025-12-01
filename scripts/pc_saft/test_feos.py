"""test_feos"""

import numpy as np


import si_units as si  # pylint: disable=E0401


from feos import Identifier, PureRecord, Parameters  # pylint: disable=E0401,E0611
from feos import EquationOfState, State, PhaseEquilibrium  # pylint: disable=E0401,E0611
from feos import PhaseDiagram  # pylint: disable=E0401,E0611


from thermolib import PcSaftGlyMix2  # pylint: disable=E0401,E0611


import matplotlib.pyplot as plt


SO2 = PureRecord(
    Identifier(name="sulfur dioxide"),
    molarweight=64.045,
    m=2.8611,
    sigma=2.6826,
    epsilon_k=205.35,
)  # sulfur_dioxide
H2O_2B = PureRecord(
    Identifier(name="water"),
    molarweight=18.015,
    m=1.0656,
    sigma=3.0007,
    epsilon_k=366.51,
    association_sites=[
        {
            "kappa_ab": 0.034868,
            "epsilon_k_ab": 2500.7,
            "na": 1,
            "nb": 1,
        }
    ],
)  # water
H2O_3B = PureRecord(
    Identifier(name="water"),
    molarweight=18.015,
    m=1.0656,
    sigma=3.0007,
    epsilon_k=366.51,
    association_sites=[
        {
            "kappa_ab": 0.034868,
            "epsilon_k_ab": 2500.7,
            "na": 1,
            "nb": 2,
        }
    ],
)  # water
H2O_4C = PureRecord(
    Identifier(name="water"),
    molarweight=18.015,
    m=1.0656,
    sigma=3.0007,
    epsilon_k=366.51,
    association_sites=[
        {
            "kappa_ab": 0.034868,
            "epsilon_k_ab": 2500.7,
            "na": 2,
            "nb": 2,
        }
    ],
)  # water
CO2_QQ = PureRecord(
    Identifier(name="carbon dioxide"),
    molarweight=44.01,
    m=1.5131,
    sigma=3.1869,
    epsilon_k=163.33,
    q=4.4,
)  # carbon_dioxide
ACETONE_DD = PureRecord(
    Identifier(name="acetone"),
    molarweight=58.08,
    m=2.7447,
    sigma=3.2742,
    epsilon_k=232.99,
    mu=2.88,
)  # acetone
CH3OH_2B = PureRecord(
    Identifier(name="methanol"),
    molarweight=32.042,
    m=1.5255,
    sigma=3.2300,
    epsilon_k=188.90,
    association_sites=[
        {
            "kappa_ab": 0.035176,
            "epsilon_k_ab": 2899.5,
            "na": 1,
            "nb": 1,
        }
    ],
)  # methanol
CH3OH_3B = PureRecord(
    Identifier(name="methanol"),
    molarweight=32.042,
    m=1.5255,
    sigma=3.2300,
    epsilon_k=188.90,
    association_sites=[
        {
            "kappa_ab": 0.035176,
            "epsilon_k_ab": 2899.5,
            "na": 1,
            "nb": 2,
        }
    ],
)  # methanol
CH4 = PureRecord(
    Identifier(name="methane"),
    molarweight=16.043,
    m=1.0000,
    sigma=3.7039,
    epsilon_k=150.03,
)  # methane (Gross and Sadowski 2001)
C2H6 = PureRecord(
    Identifier(name="ethane"),
    molarweight=30.07,
    m=1.6069,
    sigma=3.5206,
    epsilon_k=191.42,
)  # ethane (Gross and Sadowski 2001)
CO2 = PureRecord(
    Identifier(name="carbon dioxide"),
    molarweight=44.01,
    m=2.0729,
    sigma=2.7852,
    epsilon_k=169.21,
)  # carbon_dioxide
CL2_QQ = PureRecord(
    Identifier(name="chlorine"),
    molarweight=70.905,
    m=1.4682,
    sigma=3.4480,
    epsilon_k=269.67,
    q=3.0724,
)  # chlorine
C6H6_QQ = PureRecord(
    Identifier(name="benzene"),
    molarweight=78.114,
    m=2.2463,
    sigma=3.7852,
    epsilon_k=296.24,
    q=5.5907,
)  # benzene
BUTANONE_DD = PureRecord(
    Identifier(name="butanone"),
    molarweight=72.107,
    m=2.9835,
    sigma=3.4239,
    epsilon_k=244.99,
    mu=2.78,
)  # butanone
PROPANAL_DD = PureRecord(
    Identifier(name="propanal"),
    molarweight=58.08,
    m=2.6001,
    sigma=3.2872,
    epsilon_k=235.21,
    mu=2.72,
)  # propanal


def test_pc_saft_pure():  # pylint: disable=too-many-statements,disable=too-many-branches
    """test_pc_saft_pure"""
    # SO2
    parameters = Parameters.new_pure(SO2)
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
    parameters = Parameters.new_pure(H2O_2B)
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 51179:
        print("Error in tp_flash :: rho() # H20_2B")
    # H2O_3B
    parameters = Parameters.new_pure(H2O_3B)
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 51371:
        print("Error in tp_flash :: rho() # H20_3B")
    # H2O_4C
    parameters = Parameters.new_pure(H2O_4C)
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 55670:
        print("Error in tp_flash :: rho() # H20_4C")
    # CO2_QQ
    parameters = Parameters.new_pure(CO2_QQ)
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 41:
        print("Error in tp_flash :: rho() # CO2_QQ")
    # ACETONE_DD
    parameters = Parameters.new_pure(ACETONE_DD)
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 13337:
        print("Error in tp_flash :: rho() # ACEONE_DD")
    # CH3OH_2B
    parameters = Parameters.new_pure(CH3OH_2B)
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
    parameters = Parameters.new_pure(CH3OH_3B)
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
    parameters = Parameters.new_binary([CH4, C2H6])
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
    parameters = Parameters.new_binary([CH4, C2H6])
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
    parameters = Parameters.new_binary([CO2, CH3OH_2B])
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
    parameters = Parameters.new_binary([CO2_QQ, CH3OH_2B])
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
    parameters = Parameters.new_binary([ACETONE_DD, CH3OH_2B])
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
    parameters = Parameters.new_binary([CL2_QQ, C6H6_QQ])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 14137:
        print("Error in tp_flash :: rho() # Cl2_QQ+C6H6_QQ")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.97915:
        print(f"Error in tx_flash :: x = {0.5}, # Cl2_QQ+C6H6_QQ")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.02098:
        print(f"Error in ty_flash :: y = {0.5}, # Cl2_QQ+C6H6_QQ")
    # BUTANONE_DD+PROPANAL_DD
    parameters = Parameters.new_binary([BUTANONE_DD, PROPANAL_DD])
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
    parameters = Parameters.new_binary([CO2_QQ, ACETONE_DD])
    eos = EquationOfState.pcsaft(parameters)
    state = State(
        eos,
        molefracs=np.asarray([0.5, 0.5]),
        temperature=298.15 * si.KELVIN,
        pressure=0.1e6 * si.PASCAL,
    )
    if np.round(state.density / (si.MOL / si.METER**3)) != 41:
        print("Error in tp_flash :: rho() # CO2_QQ+ACETONE_DD")
    bubble = PhaseEquilibrium.bubble_point(
        eos, 298.15 * si.KELVIN, np.array([0.5, 0.5])
    )
    if np.round(bubble.vapor.molefracs[0] * 1e5) / 1e5 != 0.98973:
        print(f"Error in tx_flash :: x = {0.5}, # CO2_QQ+ACETONE_DD")
    dew = PhaseEquilibrium.dew_point(eos, 298.15 * si.KELVIN, np.array([0.5, 0.5]))
    if np.round(dew.liquid.molefracs[0] * 1e5) / 1e5 != 0.00688:
        print(f"Error in ty_flash :: y = {0.5}, # CO2_QQ+ACETONE_DD")


def test_phase_diagram():
    """test_phase_diagram"""
    fluid = PcSaftGlyMix2(
        [0.5, 0.5],
        [1.5131, 1.5255],
        [3.1869, 3.2300],
        [163.33, 188.90],
        0.0,
    )
    fluid.set_2B_assoc_term(1, 0.035176, 2899.5, 1, 1, 1)  # AssocTerm
    fluid.set_polar_term([4.4, 0.0], [4, 0])  # PolarTerm
    parameters = Parameters.new_binary([CO2_QQ, CH3OH_2B])
    eos = EquationOfState.pcsaft(parameters)  # CO2_QQ+CH3OH_2B
    temp_list = [298, 323, 348, 373, 423, 473]
    color_list = ["b", "g", "r", "c", "m", "y", "k", "w"]
    for temp, color in zip(temp_list, color_list):
        print(f"Plot temp = {temp}, # CO2_QQ+CH3OH_2B")
        pxy = np.array(fluid.t2pxy_phase_diagram(temp))
        plt.plot(pxy[:, 1], pxy[:, 0] / 1e6, color=color)
        plt.plot(pxy[:, 2], pxy[:, 0] / 1e6, color=color)
        out = PhaseDiagram.binary_vle(eos, temp * si.KELVIN)
        plt.scatter(
            np.array(out.liquid.molefracs).T[0],
            out.liquid.pressure / si.PASCAL / 1e6,
            color=color,
            s=10,
        )
        plt.scatter(
            np.array(out.vapor.molefracs).T[0],
            out.vapor.pressure / si.PASCAL / 1e6,
            color=color,
            s=10,
            label=f"{temp}",
        )
        plt.scatter(pxy[-1, 1], pxy[-1, 0] / 1e6, color=color, marker="+")
        plt.scatter(pxy[-1, 2], pxy[-1, 0] / 1e6, color=color, marker="+")
    plt.title("TestPhaseDiagram")
    plt.xlabel("X(CO2)")
    plt.ylabel("P(MPa)")
    plt.legend()
    plt.show()


def main():
    """main"""
    test_pc_saft_pure()
    test_pc_saft_mix2()
    test_phase_diagram()


if __name__ == "__main__":
    main()
