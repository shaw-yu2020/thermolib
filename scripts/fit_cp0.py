"""
Cp0
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


# CODATA2018 constants: speed of light in vacuum
C = 299792458  # m s^-1
# CODATA2018 constants: Planck constant
H = 6.62607015e-34  # J Hz^-1
# CODATA2018 constants: Boltzmann constant
K = 1.380649e-23  # J K^-1
# CODATA2018 constants: molar gas constant
R = 8.314462618  # J mol^-1 K^-1


def gaussian():
    # B3LYP/cc-pVTZ+d opt freq => sclZPE=0.9886
    # B2PLYPD3/cc-pVTZ opt freq => sclZPE=0.983
    """
    # B2PLYPD3/cc-pVTZ opt freq

    Carbon Dioxide

    0 1
     C     0.00   0.00   0.00
     O     1.16   0.00   0.00
     O    -1.16   0.00   0.00


    """
    # g16 < CO2.gjf > CO2.out
    # freq = [671.6652, 671.6652, 1371.8538, 2417.1567]
    # freq = [663.7419, 663.7419, 1344.6641, 2403.6794]
    return


def orca():
    # B3LYP cc-pVT(+d)Z opt freq => sclZPE=0.9886
    """
    ! B3LYP cc-pVT(+d)Z opt freq noautostart miniprint nopop
    %maxcore     1000
    %pal nprocs 4 end
    * xyz   0   1
    C     0.00   0.00   0.00
    O     1.16   0.00   0.00
    O    -1.16   0.00   0.00
    *
    """
    # D:\ORCA_6.0.0\orca CO2.inp > CO2.out
    # freq = [670.78, 670.78, 1370.62, 2414.32]
    return


def calc_cp0(temp, freq_plus):
    """calculate isobaric heat capacity"""
    vi = np.array([freq_plus[1:]]).transpose()
    temp = H * vi / K / temp
    expt = np.exp(-temp)
    return R * (np.sum(temp**2 * expt / (1 - expt) ** 2, axis=0) + freq_plus[0])


def aly_lee_cp0(temp, params):
    """calculate isobaric heat capacity"""
    return (
        params[0]
        + params[1] * (params[2] / temp / np.sinh(params[2] / temp)) ** 2
        + params[3] * (params[4] / temp / np.cosh(params[4] / temp)) ** 2
    )


def main():
    """main function"""
    # input parameters: nonlinear molecule or not ?
    nm = False  # CO2 is linear molecule
    # input parameters: wave length (cm^-1)
    wave = np.array([670.78, 670.78, 1370.62, 2414.32])  # CO2
    # input parameters: temperature range
    temp = list(range(200, 2000, 10))

    freq = C * wave * 100  # vibration frequencies (s^-1)(Hz)
    freq_plus = np.append([4 if nm else 3.5], freq)
    cp0 = calc_cp0(temp, freq_plus)

    result = least_squares(
        lambda params, temp, cp0: aly_lee_cp0(temp, params) - cp0,
        [1, 1, 1, 1, 1],
        args=(temp, cp0),
    )
    params = np.round(result.x * 100) / 100
    print("params:", params)
    plt.scatter(temp, cp0, marker="+", c="g")
    plt.plot(temp, aly_lee_cp0(temp, params), c="r")
    plt.show()


if __name__ == "__main__":
    main()
