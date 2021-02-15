import numpy as np
from arc import *
from scipy.constants import h, k, c

def detuning(lambda_laser):
	lambda_trans_32 = 780e-9
	lambda_trans_12 = 795e-9
	lambda_laser = 850e-9
	nu_laser = c / lambda_laser
	nu_transition_12 = c / lambda_trans_12
	nu_transition_32 = c / lambda_trans_32

	det_32 = nu_laser - nu_transition_32
	det_12 = nu_laser - nu_transition_12

	return det_12, det_32


def thermal_energy(T):
	# Return the thermal energy in MHz
	return (k * T / h) / 1e6


def light_shift_Ds(atom, laser_power, w0, det_12, det_32):
	# Find the trap depth

	rabi_D1 = atom.getRabiFrequency(
	    n1 = 5,
	    l1 = 0,
	    j1 = 0.5,
	    mj1 = 0.5,
	    n2 = 5,
	    l2 = 1,
	    j2 = 0.5,
	    q = 0,
	    laserPower = laser_power,
	    laserWaist = w0
	    )
	# print(rabi_D1 / (2*np.pi)/1e6)
	# lightshift_D1 = rabi_D1**2/4/((353-384)*1e12*2*np.pi)
	# lightshift_D1_MHz = lightshift_D1/2/np.pi*1e-6


	rabi_D2 = atom.getRabiFrequency(
	    n1 = 5,
	    l1 = 0,
	    j1 = 0.5,
	    mj1 = 0.5,
	    n2 = 5,
	    l2 = 1,
	    j2 = 1.5,
	    q = 0,
	    laserPower = laser_power,
	    laserWaist = w0
	    )

	# lightshift_D2 = rabi_D2**2/4/((353-377)*1e12*2*np.pi)
	rabi_D2 /= 2 * np.pi  # In Hz
	light_shift_P_32 = rabi_D2 ** 2 / (4 * det_32 )

	rabi_D1 /= 2 * np.pi  # In Hz
	light_shift_P_12 = rabi_D1 ** 2 / (4 * det_12)

	lightshift_D2 = 2 * light_shift_P_32 + light_shift_P_12
	lightshift_D2_MHz = lightshift_D2 / 1e6

	lightshift_D1 = 2 * light_shift_P_12 + light_shift_P_32
	lightshift_D1_MHz = lightshift_D1 / 1e6

	trap_depth = (light_shift_P_12 + light_shift_P_32) / 1e6



	# Temperature correction:
	T = 60e-6
	dE_temperature = thermal_energy(T)
	lightshift_D2 += - (det_32 / np.abs(det_32)) * dE_temperature
	lightshift_D1 += - (det_32 / np.abs(det_12)) * dE_temperature

	return lightshift_D1_MHz, lightshift_D2_MHz, trap_depth


