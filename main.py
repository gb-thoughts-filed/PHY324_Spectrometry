import Spectrometry_Calculations as sc
import scipy.optimize as spo
from scipy.stats import chi2
import matplotlib.pyplot as plt
import numpy as np

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    c = 3.0*10**8
    #Speed of Light

    h = 4.135667696*10**-15
    #Planck's Constant Electron Volts

    #==================================================================================================================

    # Mercury Wavelength Analysis

    mercury_experimental_wavelengths = [[415.2, 2.5], [445.0, 2.5], [551.1, 2.5], [581.6, 2.5]]

    mercury_experimental_wavelengths = [[x[0]*10**-9, x[1]*10**-9] for x in mercury_experimental_wavelengths]

    #print(mercury_experimental_wavelengths)

    mercury_expected_wavelengths = [[407.7837, 0.00005], [435.8328, 0.00005], [546.0735, 0.00005], [579.0663, 0.00005]]

    mercury_expected_wavelengths = [[x[0]*10**-9, x[1]*10**-9] for x in mercury_expected_wavelengths]

    #print(mercury_expected_wavelengths)

    m_exper_w_values, m_exper_uncer = sc.separate_value_and_uncertainty(mercury_experimental_wavelengths)
    m_expec_w_values, m_expec_uncer = sc.separate_value_and_uncertainty(mercury_expected_wavelengths)

    #print(m_exper_w_values)

    popt_mercury, pcov_mercury = spo.curve_fit(sc.linear, m_expec_w_values, m_exper_w_values, sigma=m_exper_uncer,
                                              absolute_sigma=True)

    print(popt_mercury)

    plt.figure(1)
    plt.errorbar(m_expec_w_values, m_exper_w_values, yerr=m_exper_uncer, marker="o", ls='',
                 label="Raw Data")
    plt.plot(m_expec_w_values, sc.linear(m_expec_w_values, *popt_mercury), marker='none', ls='-',
             label = "linear fit")
    plt.title("Theoretical Wavelengths for Mercury (m) v.s. Experimental Wavelengths for Mercury (m)")
    plt.xlabel("Theoretical Wavelengths for Mercury (m)")
    plt.ylabel("Experimental Wavelengths for Mercury (m)")
    plt.legend()

    plt.figure(2)
    plt.errorbar(m_expec_w_values, m_exper_w_values - sc.linear(m_expec_w_values, *popt_mercury), yerr=m_exper_uncer,
                 marker="o", ls='',
                 label="Difference Between Raw Data and Linear Fit")
    plt.title("Residual Plot: Theoretical Wavelengths for Mercury (m) v.s. Experimental Wavelengths for Mercury (m)")



    hg_wavelength_chi_squared = sc.reduced_chi_squared(m_exper_w_values, sc.linear(m_expec_w_values, *popt_mercury),
                                                       m_exper_uncer, len(m_exper_w_values), 2)

    hg_wavelength_chi_squared_prob = (1-chi2.cdf(hg_wavelength_chi_squared, 2))

    #print(hg_wavelength_chi_squared)
    #print(hg_wavelength_chi_squared_prob)

    # ==================================================================================================================

    # Mercury Energy Analysis

    expected_energies = sc.energies(mercury_expected_wavelengths, c, h)
    experimental_energies = sc.energies(mercury_experimental_wavelengths, c, h)

    #print(experimental_energies)
    #print(expected_energies)


    m_energy_exper_w_values, m_energy_exper_uncer = sc.separate_value_and_uncertainty(experimental_energies)
    m_energy_expec_w_values, m_energy_expec_uncer = sc.separate_value_and_uncertainty(expected_energies)

    popt_energy_mercury, pcov_energy_mercury = spo.curve_fit(sc.linear, m_energy_expec_w_values,
                                                             m_energy_exper_w_values,
                                                             sigma=m_energy_exper_uncer,
                                                             absolute_sigma=True)
    print(popt_energy_mercury)
    #print(pcov_energy_mercury)

    plt.figure(3)
    plt.errorbar(m_energy_expec_w_values, m_energy_exper_w_values, yerr=m_energy_exper_uncer, marker="o", ls='',
                 label="Data")
    plt.plot(m_energy_expec_w_values, sc.linear(m_energy_expec_w_values, *popt_energy_mercury), marker='none', ls='-',
             label="Linear Fit")
    plt.title("Theoretical Energies for Mercury (eV) v.s. Experimental Energies for Mercury (eV)")
    plt.xlabel("Theoretical Energies for Mercury (eV)")
    plt.ylabel("Experimental Energies for Mercury (eV)")
    plt.legend()

    plt.figure(4)
    plt.errorbar(m_energy_expec_w_values,
                 m_energy_exper_w_values - sc.linear(m_energy_expec_w_values, *popt_energy_mercury),
                 yerr=m_energy_exper_uncer,
                 marker="o", ls='',
                 label="Difference Between Data and Linear Fit")
    plt.title("Residual Plot: Theoretical Energies for Mercury (eV) v.s. Experimental Energies for Mercury (eV)")

    #plt.show()

    #print(experimental_energies)

    hg_energy_wavelength_chi_squared = sc.reduced_chi_squared(m_energy_exper_w_values,
                                                    sc.linear(m_energy_expec_w_values, *popt_energy_mercury),
                                                       m_energy_exper_uncer,
                                                            len(m_energy_exper_w_values), 2)

    hg_energy_wavelength_chi_squared_prob = (1 - chi2.cdf(hg_energy_wavelength_chi_squared, 2))

    #print(hg_energy_wavelength_chi_squared)
    #print(hg_energy_wavelength_chi_squared_prob)


    # =================================================================================================================

    # Hydrogen Analysis

    #qs -> hydrogen quantum state
    hqs3 = sc.hydrogen_energies_of_quantum_states(1, 3)
    hqs4 = sc.hydrogen_energies_of_quantum_states(1, 4)
    hqs5 = sc.hydrogen_energies_of_quantum_states(1, 5)

    #hee -> hydrogen expected energy
    hee3 = sc.balmer_series(1, 3)
    hee4 = sc.balmer_series(1, 4)
    hee5 = sc.balmer_series(1, 5)

    #hew -> hydrogen expected wavelength

    hew3 = (1/hee3)*c*h
    hew4 = (1/hee4)*c*h
    hew5 = (1/hee5)*c*h

    print(hqs3, hqs4, hqs5)
    print(hew3, hew4, hew5)

    hydrogen_data = [[492.9, 2.5], [658.2, 2.5]]
    hydrogen_data = [[x[0] * 10 ** -9, x[1] * 10 ** -9] for x in hydrogen_data]

    x_line_1 = [650*10**-9,670*10**-9]
    x_line_2 = [480 * 10 ** -9, 510 * 10 ** -9]
    y_line_coords = [0,0]

    plt.figure(5)

    #expected values
    plt.errorbar(hew3, 0, label="Peak 2 Theoretical", marker="o", color="green")
    plt.errorbar(hydrogen_data[1][0], 0, label="Peak 2 Experimental", marker="o", color="orange")
    plt.errorbar(hydrogen_data[1][0] + hydrogen_data[1][1], 0,
                 label="Peak 2 Experimental Plus Uncertainty",
                 marker="X", color="red")
    plt.errorbar(hydrogen_data[1][0] - hydrogen_data[1][1], 0,
                 label="Peak 2 Experimental Minus Uncertainty",
                 marker="X", color="red")
    plt.plot(x_line_1, y_line_coords, marker='', ls="-", color="black")
    plt.title("Theoretical Hydrogen Peak 2 Wavelength (m) v.s. Experimental Value (m)")
    plt.yticks([])
    plt.legend()

    plt.figure(6)
    plt.errorbar(hew4, 0, label="Peak 1 Theoretical", marker="o", color="green")
    plt.errorbar(hydrogen_data[0][0], 0, label="Peak 1 Experimental", marker="o", color="orange")
    plt.errorbar(hydrogen_data[0][0] + hydrogen_data[0][1], 0,
                 label="Peak 1 Experimental Plus Uncertainty",
                 marker="X", color="red")
    plt.errorbar(hydrogen_data[0][0] - hydrogen_data[0][1], 0,
                 label="Peak 1 Experimental Minus Uncertainty",
                 marker="X", color="red")
    plt.plot(x_line_2, y_line_coords, marker='', ls="-", color="black")
    plt.title("Theoretical Hydrogen Peak 1 Wavelength (m) v.s. Experimental Value (m)")
    plt.legend()
    plt.yticks([])

    #plt.show()

    #==================================================================================================================

    #Helium Analysis

    helium_data =[[399.5, 2.5], [455.6, 2.5], [506.9, 2.5], [589.4, 2.5], [668.9, 2.5], [706.3, 2.5], [727.2, 2.5]]

    lambda_r = [58.43339, 706.5190, 388.8648, 587.5621, 667.8151, 501.56783, 447.14802]
    lambda_r = [x * 10 ** -9 for x in lambda_r]
    lambda_r_energies = [h*c/x for x in lambda_r]

    heqs2 = sc.hydrogen_energies_of_quantum_states(2, 2)
    heqs3 = sc.hydrogen_energies_of_quantum_states(2, 3)


    #Transition expected energy between m=2 and n=3
    he_ee3 = sc.balmer_series(2, 3)

    #Transition expected wavelength between m=2 and n=3
    he_ew3 = (1/he_ee3)*c*h

    print(he_ew3)
    print(lambda_r_energies)











