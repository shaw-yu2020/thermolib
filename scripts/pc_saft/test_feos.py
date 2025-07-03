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


def main():
    """main"""
    test_pc_saft_gly_pure()
    test_pc_saft_pure()


if __name__ == "__main__":
    main()
