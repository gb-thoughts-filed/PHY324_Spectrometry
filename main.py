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

    mercury_experimental_wavelengths = [[415.2, 3], [445.0, 2], [551.1, 3], [581.6, 4]]

    mercury_experimental_wavelengths = [[x[0]*10**-9, x[1]*10**-9] for x in mercury_experimental_wavelengths]

    print("Mercury Experimental Wavelengths",mercury_experimental_wavelengths)

    mercury_expected_wavelengths = [[407.7837, 0.00005], [435.8328, 0.00005], [546.0735, 0.00005], [579.0663, 0.00005]]

    mercury_expected_wavelengths = [[x[0]*10**-9, x[1]*10**-9] for x in mercury_expected_wavelengths]

    print("Mercury Expected Wavelengths",mercury_expected_wavelengths)

    m_exper_w_values, m_exper_uncer = sc.separate_value_and_uncertainty(mercury_experimental_wavelengths)
    m_expec_w_values, m_expec_uncer = sc.separate_value_and_uncertainty(mercury_expected_wavelengths)

    #print(m_exper_w_values)

    # Curve fit for mercury wavelengths done here!

    popt_mercury, pcov_mercury = spo.curve_fit(sc.linear, m_expec_w_values, m_exper_w_values, sigma=m_exper_uncer,
                                              absolute_sigma=True)

    print("popt mercury",popt_mercury)
    print("pcov mercury", pcov_mercury)

    p_sigma_mercury = np.sqrt(np.diag(pcov_mercury))

    print("p sigma mercury", p_sigma_mercury)

    plt.figure(1)
    plt.errorbar(m_expec_w_values, m_exper_w_values, yerr=m_exper_uncer, marker="o", ls='',
                 label="Raw Data", markersize=10)
    plt.plot(m_expec_w_values, sc.linear(m_expec_w_values, *popt_mercury), marker='none', ls='-',
             label="Linear Fit")
    plt.title("Theoretical Wavelengths for Mercury (m) v.s. "
              "Experimental Wavelengths for Mercury (m)", wrap=True,
              fontsize=36)
    plt.xlabel("Theoretical Wavelengths for Mercury (m)", fontsize=32, wrap=True)
    plt.xticks(fontsize=18)
    plt.ylabel("Experimental Wavelengths for Mercury (m)", fontsize=32, wrap=True)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=24)

    plt.figure(2)
    plt.errorbar(m_expec_w_values, m_exper_w_values - sc.linear(m_expec_w_values, *popt_mercury), yerr=m_exper_uncer,
                 marker="o", ls='', markersize=10,
                 label="Difference Between Raw Data and Linear Fit")
    plt.title("Residual Plot: Theoretical Wavelengths for Mercury (m) "
              "v.s. Experimental Wavelengths for Mercury (m)",
              wrap=True, fontsize=36)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Theoretical Wavelengths for Mercury (m)",
               fontsize=32, wrap=True)
    plt.ylabel("Difference Between Observed Experimental Wavelengths"
               " and Linear Fit (m)", wrap=True, fontsize=32)
    plt.axhline(y=0, color="red", ls="--")
    plt.legend(fontsize=24, loc=8)


    hg_wavelength_chi_squared = sc.reduced_chi_squared(m_exper_w_values, sc.linear(m_expec_w_values, *popt_mercury),
                                                       m_exper_uncer, len(m_exper_w_values), 2)

    hg_wavelength_chi_squared_prob = (1-chi2.cdf(hg_wavelength_chi_squared, 2))

    print("Mercury Chi-Squared", hg_wavelength_chi_squared)
    print("Mercury Chi-Squared Probability" ,hg_wavelength_chi_squared_prob)

    # ==================================================================================================================

    # Mercury Energy Analysis

    expected_energies = sc.energies(mercury_expected_wavelengths, c, h)
    energy_one = sc.energies([[(404.6565*10**-9), (0.00005*10**-9)]], c, h)
    experimental_energies = sc.energies(mercury_experimental_wavelengths, c, h)

    print("Mercury Experimental Energies", experimental_energies)
    print("Mercury Expected Energies", energy_one, expected_energies)


    m_energy_exper_w_values, m_energy_exper_uncer = sc.separate_value_and_uncertainty(experimental_energies)
    m_energy_expec_w_values, m_energy_expec_uncer = sc.separate_value_and_uncertainty(expected_energies)

    popt_energy_mercury, pcov_energy_mercury = spo.curve_fit(sc.linear, m_energy_expec_w_values,
                                                             m_energy_exper_w_values,
                                                             sigma=m_energy_exper_uncer,
                                                             absolute_sigma=True)
    print("popt energy mercury", popt_energy_mercury)
    print("pcov energy mercury", pcov_energy_mercury)

    p_sigma_mercury_energy = np.sqrt(np.diag(pcov_energy_mercury))
    print("p sigma energy mercury", p_sigma_mercury_energy)

    plt.figure(3)
    plt.errorbar(m_energy_expec_w_values, m_energy_exper_w_values, yerr=m_energy_exper_uncer, marker="o", ls='',
                 label="Data", markersize=10)
    plt.plot(m_energy_expec_w_values, sc.linear(m_energy_expec_w_values, *popt_energy_mercury), marker='none', ls='-',
             label="Linear Fit")
    plt.title("Theoretical Energies for Mercury (eV) v.s. "
              "Experimental Energies for Mercury (eV)", wrap=True, fontsize=36)
    plt.xlabel("Theoretical Energies for Mercury (eV)",
               wrap=True, fontsize=32)
    plt.ylabel("Experimental Energies for Mercury (eV)",
               wrap=True, fontsize=32)
    plt.legend(fontsize=24)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.figure(4)
    plt.errorbar(m_energy_expec_w_values,
                 m_energy_exper_w_values - sc.linear(m_energy_expec_w_values, *popt_energy_mercury),
                 yerr=m_energy_exper_uncer,
                 marker="o", ls='', markersize=10,
                 label="Difference Between Data and Linear Fit")
    plt.title("Residual Plot: Theoretical Energies for Mercury (eV) v.s. "
              "Experimental Energies for Mercury (eV)", wrap=True,
              fontsize=36)
    plt.axhline(y=0, color="red", ls="--")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Theoretical Energies for Mercury (eV)",
               fontsize=32, wrap=True)
    plt.ylabel("Difference Between Calculated Energies Based On Observed "
               "Wavelengths"
               " and Linear Fit (eV)", wrap=True, fontsize=32)
    plt.legend(fontsize=24, loc=8)
    #plt.show()

    #print(experimental_energies)

    hg_energy_wavelength_chi_squared = sc.reduced_chi_squared(m_energy_exper_w_values,
                                                    sc.linear(m_energy_expec_w_values, *popt_energy_mercury),
                                                       m_energy_exper_uncer,
                                                            len(m_energy_exper_w_values), 2)

    hg_energy_wavelength_chi_squared_prob = (1 - chi2.cdf(hg_energy_wavelength_chi_squared, 2))

    print("Mercury Energy Chi-Squared",hg_energy_wavelength_chi_squared)
    print("Mercury Energy Chi-Squared Probability",
          hg_energy_wavelength_chi_squared_prob)


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

    hew = [hew3, hew4, hew5]


    print("Hydrogen Quantum States", hqs3, hqs4, hqs5)
    print("Hydrogen Expected Energies", hee3, hee4, hee5)
    print("Hydrogen Expected Wavelength", hew3, hew4, hew5)

    hydrogen_data = [[492.9, 6], [658.2, 6]]
    hydrogen_data = [[x[0] * 10 ** -9, x[1] * 10 ** -9] for x in hydrogen_data]
    hydrogen_data_calibrated = [(x[0] - popt_mercury[1])/popt_mercury[0] for x in hydrogen_data]
    print("Hydrogen Data Calibrated", hydrogen_data_calibrated)


    plt.figure(5)

    #expected values
    plt.plot(hew3, 0, label="Peak 2 Theoretical",
                 marker="o", color="green", markersize=22, linewidth=0)
    plt.plot(hydrogen_data_calibrated[1], 0, label="Peak 2 Experimental",
                 marker="o", color="orange", markersize=22, linewidth=0)
    plt.plot(hydrogen_data_calibrated[1] + hydrogen_data[1][1], 0,
                 label="Peak 2 Experimental Plus Uncertainty",
                 marker="X", color="red", markersize=22, linewidth=0)
    plt.plot(hydrogen_data_calibrated[1] - hydrogen_data[1][1], 0,
                 label="Peak 2 Experimental Minus Uncertainty",
                 marker="X", color="red", markersize=22, linewidth=0)
    plt.axhline(y=0, color="black", ls="-")
    plt.title("Theoretical Hydrogen Peak 2 Wavelength (m) "
              "v.s. Experimental Value (m)", fontsize=36, wrap=True)
    plt.xlabel("Wavelength (m)", fontsize=32)
    plt.legend(fontsize=24)
    plt.xticks(fontsize=22)
    plt.yticks([])

    plt.figure(6)
    plt.errorbar(hew4, 0, label="Peak 1 Theoretical",
                 marker="o", color="green", markersize=22, linewidth=0)
    plt.errorbar(hydrogen_data_calibrated[0], 0, label="Peak 1 Experimental",
                 marker="o", color="orange", markersize=22, lw=0)
    plt.errorbar(hydrogen_data_calibrated[0] + hydrogen_data[0][1], 0,
                 label="Peak 1 Experimental Plus Uncertainty",
                 marker="X", color="red", markersize=22, lw=0)
    plt.errorbar(hydrogen_data_calibrated[0] - hydrogen_data[0][1], 0,
                 label="Peak 1 Experimental Minus Uncertainty",
                 marker="X", color="red", markersize=22, lw=0)
    plt.axhline(y=0, color="black", ls="-")
    plt.title("Theoretical Hydrogen Peak 1 Wavelength (m) v.s. "
              "Experimental Value (m)", wrap=True, fontsize=36)
    plt.xlabel("Wavelength (m)", fontsize=32)
    plt.legend(fontsize=24)
    plt.xticks(fontsize=22)
    plt.yticks([])

    plt.show()

    #==================================================================================================================

    #Helium Analysis

    helium_data =[[399.5, 2.5], [455.6, 2.5], [506.9, 2.5], [589.4, 2.5], [668.9, 2.5], [706.3, 2.5], [727.2, 2.5]]
    helium_data = [[x[0] * 10 ** -9, x[1] * 10 ** -9] for x in helium_data]
    helium_data_calibrated = [(x[0] - popt_mercury[1])/popt_mercury[0] for x in helium_data]
    print("Calibrated Helium Data", helium_data_calibrated)
    lambda_r = [58.43339, 706.5190, 388.8648, 587.5621, 667.8151, 501.56783, 447.14802]
    lambda_r = [x * 10 ** -9 for x in lambda_r]
    lambda_r_energies = [h*c/x for x in lambda_r]

    heqs2 = sc.hydrogen_energies_of_quantum_states(2, 2)
    heqs3 = sc.hydrogen_energies_of_quantum_states(2, 3)


    #Transition expected energy between m=2 and n=3
    he_ee3 = sc.balmer_series(2, 3)

    #Transition expected wavelength between m=2 and n=3
    he_ew3 = (1/he_ee3)*c*h

    print("Helium quantum states",heqs2, heqs3)
    print("Helium and Hydrogen Spectral Line Wavelength", he_ew3, hew3)
    print("Lab Manual lambda_r energies", lambda_r_energies)


    #==========================================================================

    #Determining the Unknown Gas

    unknown_gas_data = [[762.1, 2.5], [771.5, 2.5],
                        [814.8, 2.5], [824.4, 2.5], [831.1, 2.5]]

    unknown_gas_data = [[x[0] * 10 ** -9, x[1] * 10 ** -9]
                        for x in unknown_gas_data]
    print("Unknown Gas Data", unknown_gas_data)

    unknown_gas_data_calibrated = [(x[0] - popt_mercury[1]*1.25)/popt_mercury[0]
                                   for x in unknown_gas_data]

    print("Unknown Gas Data Calibrated", unknown_gas_data_calibrated)

    unknown_gas_data_2 = [[561.9, 2.5], [590.8, 2.5], [761.1, 2.5],
                          [771.3, 2.5], [784.8, 2.5], [814.8, 2.5], [830.6, 2.5]]

    unknown_gas_data_2 = [[x[0] * 10 ** -9, x[1] * 10 ** -9]
                        for x in unknown_gas_data_2]


    print("Unknown Gas Data 2", unknown_gas_data_2)

    unknown_gas_data_2_calibrated = [(x[0] - popt_mercury[1]*1.25)/popt_mercury[0]
                                   for x in unknown_gas_data_2]

    print("Unknown Gas Data 2 Calibrated", unknown_gas_data_2_calibrated)

    unknown_gas_intensity_data = [4.100, 4.458, 30.775, 11.984, 6.910,
                                  26.247, 16.830]
    max_intensity = unknown_gas_intensity_data[2]
    unknown_gas_relative_intensities = [(x/max_intensity)*100 for x in
                                        unknown_gas_intensity_data]

    print("Unknown Gas Intensities", unknown_gas_relative_intensities)








