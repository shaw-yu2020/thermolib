"""giiz_gly"""

import numpy as np
import matplotlib.pyplot as plt


def gii(eta, f0, f1, f2):
    """gii"""
    return (
        f0 / (1 - eta)
        + 1.5 * f1 * eta / (1 - eta) ** 2
        + 0.5 * f2 * eta**2 / (1 - eta) ** 3
    ) / f0


def main():
    """main"""
    d = np.array(range(0, 7504)) / 10000
    g = gii(d, 1.0714, 1.2621, -0.3698) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="2B-methanol")
    g = gii(d, 1.1929, 0.6080, 0.5116) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="3B-methanol")
    g = gii(d, 1.0025, 0.9962, 0.9970) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="2B-ethanol")
    g = gii(d, 1.3296, 0.4443, 1.4808) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="3B-ethanol")
    g = gii(d, 1.4482, 0.1625, 2.2331) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="2B-npropanol")
    g = gii(d, 0.4834, 0.0918, 3.9092) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="3B-npropanol")
    g = gii(d, 1.1767, -0.2650, 4.4788) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="2B-ipropanol")
    g = gii(d, 1.2635, -0.1568, 3.4404) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="3B-ipropanol")
    g = gii(d, 1.8210, -1.2714, 6.5992) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="2B-nbutanol")
    g = gii(d, 1.8784, -0.7116, 4.8629) / gii(d, 1, 1, 1)
    plt.plot(d, g, label="3B-nbutanol")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
