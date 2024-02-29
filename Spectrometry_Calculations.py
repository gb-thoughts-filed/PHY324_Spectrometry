import PHY324_Module
import numpy as np
def linear(x, a, b):

    return a*x + b

def non_linear(x, c, d, e):

    return c+x++d +e

def energies(expected_wavelength_lst: list, planck_const, speed_of_light):

    energies = []

    for i in expected_wavelength_lst:
        energies.append(

            [planck_const*speed_of_light/i[0], planck_const*speed_of_light/i[1]]

        )

    return energies

def separate_value_and_uncertainty(list_of_lists: list):
    lst1 = []
    lst2 = []

    for i in list_of_lists:
        lst1.append(i[0])
        lst2.append(i[1])

    return np.array(lst1), np.array(lst2)
