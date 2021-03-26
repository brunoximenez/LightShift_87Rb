import numpy as np
from arc import *
from scipy.constants import h, hbar, k, c, epsilon_0, e
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize, fit_report
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


def build_params():
	params = Parameters()
	params.add_many(
		('a0', 1.),
		('a1', 1.),
		('a2', 0.1),
		)
	return params


def model_pol(params, x, data=[]):
	a0 = params['a0'].value
	a1 = params['a1'].value
	a2 = params['a2'].value

	model = a0 * x + a1 * x ** 2 + a2 * x ** 3

	if len(data) == 0:
		return model
	else:
		return model - data


data_cal_v_power = np.loadtxt('VtoP_data.dat')
plt.plot(data_cal_v_power[:, 0], data_cal_v_power[:, 1], 'o', markerfacecolor='blue')

params = build_params()
fit_output = minimize(model_pol, params, args=(data_cal_v_power[:, 0], data_cal_v_power[:, 1]), method='leaastsq')

bestfit_x = np.linspace(0, 10, 50)
bestfit_y = model_pol(fit_output.params, bestfit_x)
plt.plot(bestfit_x, bestfit_y, 'r--')
print(fit_report(fit_output.params))
fontsize = 14
plt.xlabel('Sequence channel output (V)', fontsize=fontsize)
plt.ylabel('Trap power (mW)', fontsize=fontsize)
plt.ylim((0, 5))
plt.xlim((0, np.amax(data_cal_v_power[:, 0])))


def model_linear(params, x, data=[]):
	a0 = params['a0'].value
	a1 = params['a1'].value

	model = a1 * x + a0

	if len(data) == 0:
		return model
	else:
		return model - data


def build_params_linear():
	params = Parameters()
	params.add_many(
		('a0', 1.),
		('a1', 1.),
		)
	return params


data_cal_v_ls = np.loadtxt('VtoLS_data.dat')
calibrated_trap_power = model_pol(fit_output.params, data_cal_v_ls[:, 0])
plt.figure()


params_lightshift = build_params_linear()
fit_output_lightshift = minimize(model_linear, params_lightshift, args=(calibrated_trap_power, data_cal_v_ls[:, 1]), method='leastsq')
print(fit_report(fit_output_lightshift.params))


best_fit_x_ls = np.linspace(0, 5, 50)
best_fit_y_ls = model_linear(fit_output_lightshift.params, best_fit_x_ls)

plt.plot(best_fit_x_ls, best_fit_y_ls - fit_output_lightshift.params['a0'].value, 'r--', label='Best fit')
plt.plot(calibrated_trap_power, data_cal_v_ls[:, 1] - fit_output_lightshift.params['a0'].value, 'o', markerfacecolor='blue', label='Data 6P j=3/2 mf=3')
plt.ylim((0, 25))
plt.xlim((0, 5))
plt.xlabel('Trap power (mW)', fontsize=fontsize)
plt.ylabel('Light shift (MHz)', fontsize=fontsize)



'''
Calculating now the theoretical value of light shift for the 6P transition
'''
# laser_power = 4.e-3
laser_wavelength = 850e-9
w0 = 1.1e-6
q = -1

n = 5
l = 0
j = .5
f = 2
mf = 2
state_to_calculate_shift = np.array([n, l, j, f, mf])

# Which manifold to consider? n, l, j
coupled_states = np.array([[5, 1, .5], [5, 1, 1.5]])
laser_power = np.linspace(0, 5, 6)
ls_ground_vector = np.array([])

for trap_power in laser_power:
	# Now loop over them and calculate light shift on state to calculate for each of these states to consider
	ls_ground = 0 
	for states in coupled_states:
	    # Get complete HF manifold
	    hf_manifold = ls.states_assemble(states[0], states[1], states[2])
	    # Calculate shift fue to hf_manifold
	    ls_ground += ls.light_shift_calc(state_to_calculate_shift,
	                      hf_manifold, q, trap_power * 1e-3, w0, laser_wavelength)

	# print('Light on ground state: ', ls_ground)
	ls_ground_vector = np.append(ls_ground_vector, ls_ground)
print(ls_ground_vector)
'''
Now for the excited state
'''
n = 6
l = 1
j = 1.5
f = 3
mf = 3
state_to_calculate_shift = np.array([n, l, j, f, mf])

# Which manifold to consider? n, l, j
coupled_states = np.array([[6, 0, .5], [7, 0, .5], [8, 0, .5], [5, 2, 1.5]])
ls_excited_vector = np.array([])

for trap_power in laser_power:
	# Now loop over them and calculate light shift on state to calculate for each of these states to consider
	ls_excited = 0 
	for states in coupled_states:
	    # Get complete HF manifold
	    hf_manifold = ls.states_assemble(states[0], states[1], states[2])
	    # Calculate shift fue to hf_manifold
	    ls_excited += ls.light_shift_calc(state_to_calculate_shift,
	                      hf_manifold, q, trap_power * 1e-3, w0, laser_wavelength)
	    print('Light shift due to state (n,l,j): ', states)
	    print(ls.light_shift_calc(state_to_calculate_shift,
	                      hf_manifold, q, trap_power * 1e-3, w0, laser_wavelength))

	# print('Light on excited state: ', ls_excited)
	ls_excited_vector = np.append(ls_excited_vector, ls_excited)
print(ls_excited_vector)

ls_transition = ls_excited_vector - ls_ground_vector
print('ls on transition', ls_transition)

plt.plot(laser_power, ls_transition, '--', color='black', label='Theoretical value')
plt.legend()
plt.show()


