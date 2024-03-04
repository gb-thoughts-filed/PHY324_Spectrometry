import numpy as np
from matplotlib import pyplot as plt

from PHY324_Module import *


def calculate_absorbance_transmittance(T_measured, T_uncertainty):
    A_measured = -np.log10(T_measured)

    A_uncertainty = np.abs(-1 / (T_measured * np.log(10))) * T_uncertainty

    return A_measured, A_uncertainty


# Example usage:
if __name__ == "__main__":
    T_measured = 0.458 / 100
    T_uncertainty = 0.000005  # Uncertainty in transmittance

    results = round_to_error(*calculate_absorbance_transmittance(T_measured, T_uncertainty))
    print(results)
    plt.figure("Residual for A")
    t_colour_table = {"blue": [[63.23, 2.52], [0.199, 1.598]], "red": [[96.485, 0.458], [0.015, 2.6]],
                      "green": [[36.34, 7.778], [0.44, 1.109]]}
    for colour in t_colour_table:
        a_calc = []
        for i in t_colour_table[colour][0]:
            a_calc.append(calculate_absorbance_transmittance(i / 100, 0.000005)[0])
        plt.errorbar(t_colour_table[colour][0], np.subtract(t_colour_table[colour][1], a_calc),
                     yerr=np.full_like(t_colour_table[colour][0], 0.0005),
                     label='Residual for Absorbance of {} coloured dye'.format(colour),
                     ls='',
                     lw=1,
                     marker='o', mfc=colour, mec=colour,
                     markersize=5, capsize=1.5, capthick=0.5, zorder=0, ecolor="grey", alpha=0.9)
    plt.axhline(y=0, color='red', linestyle='--', label='Zero Error')
    plt.legend()
    plt.title("Residuals of Absorbance for Different Transmittance")
    plt.xlabel("Measured Transmittance (%)")
    plt.ylabel("Residuals of Absorbance (Measured - Calculated)")
    plt.tight_layout()
    plt.savefig("residual_A_full")

    plt.figure("Residual for A with out outlier")
    t_colour_table = {"blue": [[63.23, 2.52], [0.199, 1.598]], "red": [[96.485], [0.015]],
                      "green": [[36.34, 7.778], [0.44, 1.109]]}
    for colour in t_colour_table:
        a_calc = []
        for i in t_colour_table[colour][0]:
            a_calc.append(calculate_absorbance_transmittance(i / 100, 0.000005)[0])
        plt.errorbar(t_colour_table[colour][0], np.subtract(t_colour_table[colour][1], a_calc),
                     yerr=np.full_like(t_colour_table[colour][0], 0.0005),
                     label='Residual for Absorbance of {} coloured dye'.format(colour),
                     ls='',
                     lw=1,
                     marker='o', mfc=colour, mec=colour,
                     markersize=5, capsize=1.5, capthick=0.5, zorder=0, ecolor="grey", alpha=0.9)
    plt.axhline(y=0, color='red', linestyle='--', label='Zero Error')
    plt.legend()
    plt.title("Residuals of Absorbance for Different Transmittance\nexcluding T=0.458")
    plt.xlabel("Measured Transmittance (%)")
    plt.ylabel("Residuals of Absorbance (Measured - Calculated)")
    plt.tight_layout()
    plt.savefig("residual_A_partial")
    plt.show()
