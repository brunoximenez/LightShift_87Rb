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

a0 = 5.29177210903e-11
fontsize = 16
atom = Rubidium87()

def light_shift_calc(state_to_calculate_shift, cooupled_states, q, laser_power, w0, laser_wavelength):
    light_shift = 0
    laser_frequency = c / laser_wavelength

    n = int(state_to_calculate_shift[[0]])
    l = int(state_to_calculate_shift[[1]])
    j = float(state_to_calculate_shift[[2]])
    f = float(state_to_calculate_shift[[3]])
    mf = float(state_to_calculate_shift[[4]])

    # print('\nCalculating light shift on state: ', n, l, j, f, mf)

    for n_prime, l_prime, j_prime, f_prime, mf_prime in cooupled_states:
        n_prime = int(n_prime)
        l_prime = int(l_prime)


        # if mf_prime == mf + q and abs(f_prime - f) < 2:
        if abs(f_prime - f) < 2:
            transition_frequency = atom.getTransitionFrequency(
                n, l, j, n_prime, l_prime, j_prime, s=0.5, s2=None)
            
            # print('\nCoupled state: ', n_prime, l_prime, j_prime, f_prime, mf_prime)
            # print('Wavelength of the transition (nm): ', c * 1e9 / transition_frequency)
            # print('Detuning (THz): ', (np.abs(transition_frequency) - laser_frequency) * 1e-12)
            q=1
            dip_elem = np.sqrt(.5) * atom.getDipoleMatrixElementHFS(
                n, l, j, f, mf, n_prime, l_prime, j_prime, f_prime, mf_prime, 1) * a0 * e
            dip_elem += np.sqrt(.5) * atom.getDipoleMatrixElementHFS(
                n, l, j, f, mf, n_prime, l_prime, j_prime, f_prime, mf_prime, -1) * a0 * e

            dip_elem_square = np.square(np.abs(dip_elem)) 
            # dip_elem +=  atom.getDipoleMatrixElementHFS(n, l, j, f, mf, n_prime, l_prime, j_prime, f_prime, mf_prime, -1) * \
            # 	atom.getDipoleMatrixElementHFS(n_prime, l_prime, j_prime, f_prime, mf_prime, n, l, j, f, mf, 1)

            # dip_elem +=  atom.getDipoleMatrixElementHFS(n, l, j, f, mf, n_prime, l_prime, j_prime, f_prime, mf_prime, 1) * \
            # 	atom.getDipoleMatrixElementHFS(n_prime, l_prime, j_prime, f_prime, mf_prime, n, l, j, f, mf, -1)

            # print('Dipole matrix element (au): ', dip_elem)

            # dip_elem_square = dip_elem * ((a0 * e) ** 2)
            # print('dipole matrix elemen: ', np.sqrt(dip_elem_square))


            light_shift_tmp = - dip_elem_square * \
                ls.get_E_square(laser_power, w0) / (4 * h)
            light_shift_tmp *=  (1 / (transition_frequency - laser_frequency) +
                                1 / (transition_frequency + laser_frequency))

            # In MHz
            light_shift_tmp /= h * 1e6

            print('Shift (MHz): ', light_shift_tmp)

            light_shift += light_shift_tmp

    return light_shift


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
plt.plot(data_cal_v_power[:, 0], data_cal_v_power[:, 1], 'o', markerfacecolor='skyblue', markeredgecolor='blue')

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

# plt.plot(best_fit_x_ls, best_fit_y_ls - fit_output_lightshift.params['a0'].value, 'r--', label='Best fit')
plt.plot(calibrated_trap_power, data_cal_v_ls[:, 1] - fit_output_lightshift.params['a0'].value, 'o', markerfacecolor='skyblue', markeredgecolor='blue', label='Data 6P j=3/2 mf=3')
# plt.ylim((0, 25))
plt.xlim((0, 5))
plt.xlabel('Trap power (mW)', fontsize=fontsize)
plt.ylabel('Light shift (MHz)', fontsize=fontsize)



'''
Calculating now the theoretical value of light shift for the 6P transition
'''
# laser_power = 4.e-3
laser_wavelength = 850e-9
w0 = 1.19e-6
q = 1

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
	    ls_ground += light_shift_calc(state_to_calculate_shift,
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

coupled_states = np.array([[6, 0, .5], [7, 0, .5], [8, 0, .5], [5, 2, 1.5], [5, 2, 2.5]])
ls_excited_vector = np.array([])


for trap_power in laser_power:
	# Now loop over them and calculate light shift on state to calculate for each of these states to consider
	ls_excited = 0 
	for states in coupled_states:
	    # Get complete HF manifold
	    hf_manifold = ls.states_assemble(states[0], states[1], states[2])
	    # Calculate shift fue to hf_manifold
	    ls_excited += light_shift_calc(state_to_calculate_shift,
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
plt.plot(laser_power, ls_ground_vector, '--', label='Ground state')
plt.plot(laser_power, ls_excited_vector, '--', label='Excited state')
plt.legend()
plt.show()


