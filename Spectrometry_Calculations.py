import PHY324_Module
import numpy as np
def linear(x, a, b):

    return a*x + b

def non_linear(x, c, d, e):

    return c*x**d +e

def energies(expected_wavelength_lst: list, planck_const, speed_of_light):

    energies = []

    for i in expected_wavelength_lst:

        power_uncertainty = np.abs((-1)*(i[0]**-2)*i[1])
        energies.append(

            [planck_const*speed_of_light/i[0], planck_const*speed_of_light*power_uncertainty]

        )

    return energies

def separate_value_and_uncertainty(list_of_lists: list):
    lst1 = []
    lst2 = []

    for i in list_of_lists:
        lst1.append(i[0])
        lst2.append(i[1])

    return np.array(lst1), np.array(lst2)

def reduced_chi_squared(data_ys, predicted_ys, uncertainties, num_elements, num_parameters):

    dof = num_elements - num_parameters
    residuals = np.square(data_ys - predicted_ys)
    squared_uncertainties = np.square(uncertainties)
    red_chi_squared = (1/dof)*(np.sum(residuals/squared_uncertainties))

    return red_chi_squared

def hydrogen_energies_of_quantum_states(Z, n):
    return -(13.6*(Z**2))/(n**2)

def balmer_series(Z, n):
    rydberg_eh = 13.605693

    return rydberg_eh*(Z**2)*((1/4)-(1/n**2))

