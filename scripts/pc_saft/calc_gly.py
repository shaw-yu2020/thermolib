"""PcSaftGlyPure"""

import math
from thermolib import PcSaftGlyPure  # pylint: disable=no-name-in-module


def g2b_methanol():
    """g2b_methanol"""
    m, sigma, epsilon = 2.1049, 2.9008, 195.80  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.06085, 2477.7  # kappa_ab,epsilon_ab
    f0, f1, f2 = 1.0714, 1.2621, -0.3698
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, f0, f1, f2)
    return fluid


def g2b_ethanol():
    """g2b_ethanol"""
    m, sigma, epsilon = 2.2266, 3.2625, 214.44  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.02457, 2601.7  # kappa_ab,epsilon_ab
    f0, f1, f2 = 1.0025, 0.9962, 0.9970
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, f0, f1, f2)
    return fluid


def g2b_npropanol():
    """g2b_npropanol"""
    m, sigma, epsilon = 1.8220, 3.8639, 276.76  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.00539, 2836.6  # kappa_ab,epsilon_ab
    f0, f1, f2 = 1.4482, 0.1625, 2.2331
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, f0, f1, f2)
    return fluid


def g2b_ipropanol():
    """g2b_ipropanol"""
    m, sigma, epsilon = 2.3181, 3.5750, 247.39  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.00904, 2450.0  # kappa_ab,epsilon_ab
    f0, f1, f2 = 1.1767, -0.2650, 4.4788
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, f0, f1, f2)
    return fluid


def g2b_nbutanol():
    """g2b_nbutanol"""
    m, sigma, epsilon = 2.4655, 3.7411, 266.14  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.00896, 2605.9  # kappa_ab,epsilon_ab
    f0, f1, f2 = 1.8210, -1.2714, 6.5992
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, f0, f1, f2)
    return fluid


def gross_methanol():
    """gross_methanol"""
    m, sigma, epsilon = 1.5255, 3.2300, 188.90  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.035176, 2899.5  # kappa_ab,epsilon_ab
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, 1, 1, 1)
    return fluid


def gross_ethanol():
    """gross_ethanol"""
    m, sigma, epsilon = 2.3827, 3.1771, 198.24  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.032384, 2653.4  # kappa_ab,epsilon_ab
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, 1, 1, 1)
    return fluid


def gross_npropanol():
    """gross_npropanol"""
    m, sigma, epsilon = 2.9997, 3.2522, 233.40  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.015268, 2276.8  # kappa_ab,epsilon_ab
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, 1, 1, 1)
    return fluid


def gross_ipropanol():
    """gross_ipropanol"""
    m, sigma, epsilon = 3.0929, 3.2085, 208.42  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.024675, 2253.9  # kappa_ab,epsilon_ab
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, 1, 1, 1)
    return fluid


def gross_nbutanol():
    """gross_nbutanol"""
    m, sigma, epsilon = 2.7515, 3.6139, 259.59  # m,sigma,epsilon
    kappa_ab, epsilon_ab = 0.006692, 2544.6  # kappa_ab,epsilon_ab
    fluid = PcSaftGlyPure(m, sigma, epsilon)
    fluid.set_2B_assoc_term(kappa_ab, epsilon_ab, 1, 1, 1)
    return fluid


def main():
    """main function"""
    fluid = g2b_methanol()
    fluid.c_flash()
    crit_t, crit_p, crit_rho = fluid.T(), fluid.p(), fluid.rho()
    for temp in range(math.floor(0.6 * crit_t), math.ceil(crit_t)):
        fluid.t_flash(temp)
        print(fluid.T_s(), fluid.p_s(), fluid.rho_v(), fluid.rho_l())
    print(crit_t, crit_p, crit_rho, crit_rho)


if __name__ == "__main__":
    main()
