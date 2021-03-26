import numpy as np
from arc import *
from scipy.constants import h, hbar, k, c, epsilon_0, e
import matplotlib.pyplot as plt


a0 = 5.29177210903e-11
fontsize = 16
atom = Rubidium87()
'''
Some plot settings
'''
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

# TeX settings
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}']  # for \text command


def light_shift_calc(state_to_calculate_shift, cooupled_states, q, laser_power, w0, laser_wavelength):
    light_shift = 0
    laser_frequency = c / laser_wavelength

    n = int(state_to_calculate_shift[[0]])
    l = int(state_to_calculate_shift[[1]])
    j = float(state_to_calculate_shift[[2]])
    f = float(state_to_calculate_shift[[3]])
    mf = float(state_to_calculate_shift[[4]])

    q = int(q)

    for n_prime, l_prime, j_prime, f_prime, mf_prime in cooupled_states:
        n_prime = int(n_prime)
        l_prime = int(l_prime)
        transition_frequency = atom.getTransitionFrequency(
            n, l, j, n_prime, l_prime, j_prime, s=0.5, s2=None)
        print('transition wavelength (nm):', c * 1e9 / transition_frequency)
        print('transition frequecy (THz):' , transition_frequency * 1e-12)


        print('detuning (THz): ', (laser_frequency - np.abs(transition_frequency)) * 1e-12)
        print('detuning 780-850: ', c * 1e-12 * (1 / (850  * 1e-9) - 1 / (780 * 1e-9)))

        if mf_prime == mf + q and abs(f_prime - f) < 2:
            print('f_prime', f_prime)
            print('mf_prime', mf_prime)
            dip_elem = atom.getDipoleMatrixElementHFS(
                n, l, j, f, mf, n_prime, l_prime, j_prime, f_prime, mf_prime, q)
            # print(dip_elem)
            dip_elem_square = (dip_elem * a0 * e) ** 2
            print('dipole matrix elemen **2: ', dip_elem_square)

        else:
            dip_elem_square = 0

        light_shift_tmp = dip_elem_square * \
            get_E_square(laser_power, w0) / (4 * h)
        light_shift_tmp *= (1 / (laser_frequency - transition_frequency) +
                            1 / (transition_frequency + laser_frequency))

        # In MHz
        light_shift_tmp /= h * 1e6

        light_shift += light_shift_tmp

    return light_shift


def get_E_square(P, w0):
    E2 = 2 * get_I(P, w0) / (epsilon_0 * c)
    return E2


def get_I(P, w0):
    I = 2 * P / (np.pi * w0**2)
    return I


def states_assemble(n, l, j):
    matrix_states = np.array([[0, 0, 0, 0, 0]])
    for k in range(int(10 * (np.abs(j - atom.I))), int(10 * (j + atom.I + 1)), 10):
        f = k / 10
        for i in range(-k, k + 10, 10):
            mf = i / 10
            state = np.array([[n, l, j, f, mf]])
            matrix_states = np.concatenate((matrix_states, state), axis=0)
    matrix_states = np.delete(matrix_states, 0, 0)

    return matrix_states


def main():

    laser_power = 1.e-3
    laser_wavelength = 850e-9
    w0 = 1.1e-6
    q = -1

    n = 6
    l = 1
    j = 1.5
    f = 3
    mf = 3
    state_to_calculate_shift = np.array([n, l, j, f, mf])

    # Which manifold to consider?
    coupled_states = states_assemble(6, 0, .5)

    ls_6P = light_shift_calc(state_to_calculate_shift,
                          coupled_states, q, laser_power, w0, laser_wavelength)
    print(ls_6P)

    n = 5
    l = 0
    j = .5
    f = 2
    mf = 2
    state_to_calculate_shift = np.array([n, l, j, f, mf])

    # Which manifold to consider?
    coupled_states = states_assemble(5, 1, .5)

    ls_6P = light_shift_calc(state_to_calculate_shift,
                          coupled_states, q, laser_power, w0, laser_wavelength)


main()
