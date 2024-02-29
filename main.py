import Spectrometry_Calculations as sc
import scipy.optimize as sp
import matplotlib.pyplot as plt

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    c = 3.0*10**8

    h = 4.135667696*10**-15

    mercury_experimental_wavelengths = [[415.2, 2.5], [445.0, 2.5], [551.1, 2.5], [581.6, 2.5]]

    mercury_expected_wavelengths = [[407.7837, 0.00005], [435.8328, 0.00005], [546.0735, 0.00005], [579.0663, 0.00005]]

    expected_energies = sc.energies(mercury_expected_wavelengths, c, h)
    experimental_energies = sc.energies(mercury_experimental_wavelengths, c, h)

    #print(expected_energies)
    #print(experimental_energies)

    m_exper_w_values, m_exper_uncer = sc.separate_value_and_uncertainty(mercury_experimental_wavelengths)
    m_expec_w_values, m_expec_uncer = sc.separate_value_and_uncertainty(mercury_expected_wavelengths)

    #print(m_exper_w_values)

    popt_mercury, pcov_mercury = sp.curve_fit(sc.linear, m_expec_w_values, m_exper_w_values, sigma=m_expec_uncer,
                                              absolute_sigma=True)

    print(popt_mercury)

    plt.figure(1)
    plt.errorbar(m_expec_w_values, m_exper_w_values, yerr=m_exper_uncer, marker="o", ls='',
                 label="Raw Data")
    plt.plot(m_expec_w_values, sc.linear(m_expec_w_values, *popt_mercury), marker='none', ls='-',
             label = "linear fit")

    plt.legend()
    plt.show()



