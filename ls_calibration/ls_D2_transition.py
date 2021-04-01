import numpy as np
from arc import *
from scipy.constants import h, hbar, k, c, epsilon_0, e
import matplotlib.pyplot as plt
import ls_calculator as ls


'''
Some plot settings
'''
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

# TeX settings
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}']  # for \text command


# trap parameters
laser_wavelength = 850e-9
laser_power = 4e-3
w0 = 1.1e-6
q = 1

def ls_calc_state(state, coupled_states, trap_power=4., q=0):
	# Now loop over them and calculate light shift on state to calculate for each of these states to consider
	light_shift = 0
	for states in coupled_states:
	    # Get complete HF manifold
	    hf_manifold = ls.states_assemble(states[0], states[1], states[2])
	    # Calculate shift fue to hf_manifold
	    light_shift += ls.light_shift_calc(state,
	                      hf_manifold, q, trap_power, w0, laser_wavelength)
	    print('Light shift due to state (n,l,j): ', states)

	
	return light_shift


def ls_calc_ground(f=2, mf=2, trap_power=4e-3, q=0):
	# For the ground state
	n = 5
	l = 0
	j = .5

	state_to_calculate_shift = np.array([n, l, j, f, mf])

	# Which manifold to consider? n, l, j
	coupled_states = np.array([[5, 1, .5], [5, 1, 1.5]])

	ls_ground = ls_calc_state(state_to_calculate_shift, coupled_states, trap_power=trap_power/2, q=q)
	# ls_ground += ls_calc_state(state_to_calculate_shift, coupled_states, trap_power=laser_power/2, q=-1) 

	return ls_ground


def ls_calc_5P_j32(f=3, mf=3, trap_power=4e-3, q=0):
	# For the excited state
	n = 5
	l = 1
	j = 1.5
	f = 3
	mf = 3
	state_to_calculate_shift = np.array([n, l, j, f, mf])

	# Which manifold to consider? n, l, j
	coupled_states = np.array([[5, 0, .5], [7, 0, .5]])

	ls_excited = ls_calc_state(state_to_calculate_shift, coupled_states, trap_power=trap_power, q=q)
	# ls_excited += ls_calc_state(state_to_calculate_shift, coupled_states, trap_power=laser_power/2, q=-1)

	return ls_excited


ls_5S = ls_calc_ground(q=1)
ls_5S += ls_calc_ground(q=-1)

print('Light on ground state: ', ls_5S)

ls_5P = ls_calc_5P_j32(q=1)
ls_5P += ls_calc_5P_j32(q=-1)

print('Light on excited state: ', ls_5P)



