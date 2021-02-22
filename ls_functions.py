import numpy as np
from arc import *
from scipy.constants import h, hbar, k, c, epsilon_0, e
import matplotlib.pyplot as plt

a0 = 5.29177210903e-11
fontsize = 16
'''
Some plot settings
'''
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

# TeX settings
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}']  # for \text command


def get_dipolematrix_hfs(atom, n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q):

	elem = atom.getDipoleMatrixElementHFS(n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q)
	return elem ** 2


def s_states(atom, f, mf, q, P, w0):
	'''
	Sum for all contributions from P 1/2
	'''
	light_shift = 0
	sum_dip_square = 0
	det_12, det_32 = detuning()
	
	n1 = 5
	l1 = 0
	j1 = 1/2
	f1 = f
	mf1 = mf

	n2 = 5
	l2 = 1
	j2 = 1/2
	f_min = int(atom.I - j2)
	f_max = int(atom.I + j2)

	print('Calculating light shift for P1/2 contribution...')
	for f2 in range(f_min, f_max + 1, 1):
		for mf2 in range(-f2, f2 + 1, 1):
			print('f1 = ', f2)
			if mf2 == mf1 + q and abs(f2 - f1) < 2:
				print('mf2 = ', mf2)
				dip_elem = atom.getDipoleMatrixElementHFS(n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q)
				print(dip_elem)
				dip_elem_square = (dip_elem * a0 * e) ** 2
				sum_dip_square += dip_elem_square
			print('New f...')


	# light_shift = sum_dip_square * get_E_square(P, w0) / (2 * h * det_12)
	# # In MHz
	# light_shift /= h * 1e6
	# print('Light Shift due to P1/2 in MHz: \n', light_shift)
	light_shift_transition = sum_dip_square * get_E_square(P, w0) / (4 * h * det_12)
	# In MHz
	light_shift_transition /= h * 1e6
	light_shift += light_shift_transition
	print('Light Shift due to P1/2 in MHz: \n', light_shift_transition)
	sum_dip_square = 0


	'''
	Now for P3/2
	'''
	print('Calculating light shift for P3/2 contribution...')
	j2 = 3/2
	f_min = int(atom.I - j2)
	f_max = int(atom.I + j2)

	for f2 in range(f_min, f_max + 1, 1):
		print('f2 = ', f2)
		for mf2 in range(-f2, f2 + 1, 1):
			if mf2 == mf1 + q and abs(f1 - f2) < 2:
				print('mf2 = ', mf2)
				dip_elem = atom.getDipoleMatrixElementHFS(n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q)
				dip_elem_square = (dip_elem * a0 * e) ** 2
				sum_dip_square += dip_elem_square


	# light_shift = sum_dip_square * get_E_square(P, w0) / (2 * h * det_32)
	# # In MHz
	# light_shift /= h * 1e6
	# print('Light Shift due to P1/2 + P3/2 in MHz: \n', light_shift)

	light_shift_transition = sum_dip_square * get_E_square(P, w0) / (4 * h * det_32)
	# In MHz
	light_shift_transition /= h * 1e6
	light_shift += light_shift_transition
	print('Adding Light Shift due to P3/2 in MHz: \n', light_shift_transition)


	return light_shift


def p_3_2_states(atom, f, mf, q, P, w0):
	'''
	Extra Transitions to consider:
	---> 5D3/2 ; 775.94 nm
	---> 5D5/2 ; 775.77 nm
	---> 7S1/2 ; 740.82 nm
	'''

	'''
	Sum for all contributions from S 1/2
	'''
	sum_dip_square = 0
	light_shift = 0
	det_D12, det_D32 = detuning()
	
	n1 = 5
	l1 = 0
	j1 = 1/2


	n2 = 5
	l2 = 1
	j2 = 3/2
	f2 = f
	mf2 = mf

	f_min = int(atom.I - j1)
	f_max = int(atom.I + j1)

	print('\nCalculating light shift for 5S1/2 contribution...')
	for f1 in range(f_min, f_max + 1, 1):
		for mf1 in range(-f1, f1 + 1, 1):
			# print('f1 = ', f1)
			if mf2 == mf1 + q and abs(f2 - f1) < 2:
				# print('mf1 = ', mf1)
				dip_elem = atom.getDipoleMatrixElementHFS(n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q)
				# print(dip_elem)
				dip_elem_square = (dip_elem * a0 * e) ** 2
				sum_dip_square +=  - dip_elem_square
			# print('New f...')


	light_shift_transition = sum_dip_square * get_E_square(P, w0) / (4 * h * det_D32)
	# In MHz
	light_shift_transition /= h * 1e6
	light_shift += light_shift_transition
	print('Light Shift due to S1/2 in MHz: \n', light_shift_transition)
	sum_dip_square = 0

	'''
	Sum for all contributions from 5P3/2 ---> 5D3/2 ; 775.94 nm
	'''

	lambda_trans_5D32 = 775.94e-9
	lambda_laser = 850e-9
	nu_laser = c / lambda_laser
	nu_transition_5D32 = c / lambda_trans_5D32
	det_5D32 = nu_laser - nu_transition_5D32
	
	n1 = 5
	l1 = 1
	j1 = 3/2
	f1 = f
	mf1 = mf

	n2 = 5
	l2 = 2
	j2 = 3/2

	f_min = int(atom.I - j2)
	f_max = int(atom.I + j2)

	print('\nCalculating light shift for 5D3/2 contribution...')
	for f2 in range(f_min, f_max + 1, 1):
		for mf2 in range(-f2, f2 + 1, 1):
			# print('f1 = ', f1)
			if mf2 == mf1 + q and abs(f2 - f1) < 2:
				print('mf1 = ', mf1)
				dip_elem = atom.getDipoleMatrixElementHFS(n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q)
				print(dip_elem)
				dip_elem_square = (dip_elem * a0 * e) ** 2
				sum_dip_square += dip_elem_square
			# print('New f...')


	light_shift_transition = sum_dip_square * get_E_square(P, w0) / (4 * h * det_D32)
	# In MHz
	light_shift_transition /= h * 1e6
	light_shift += light_shift_transition
	print('Light Shift due to 5D3/2 in MHz: \n', light_shift_transition)
	sum_dip_square = 0

	'''
	Sum for all contributions from 5P3/2 ---> 5D5/2 ; 775.77 nm
	'''

	lambda_trans_5D52 = 775.77e-9
	lambda_laser = 850e-9
	nu_laser = c / lambda_laser
	nu_transition_5D52 = c / lambda_trans_5D52
	det_5D52 = nu_laser - nu_transition_5D52



	n2 = 5
	l2 = 2
	j2 = 5/2

	f_min = int(atom.I - j2)
	f_max = int(atom.I + j2)

	print('\nCalculating light shift for 5P3/2 ---> 5D5/2 contribution...')
	for f2 in range(f_min, f_max + 1, 1):
		for mf2 in range(-f2, f2 + 1, 1):
			# print('f1 = ', f1)
			if mf2 == mf1 + q and abs(f2 - f1) < 2:
				# print('mf1 = ', mf1)
				dip_elem = atom.getDipoleMatrixElementHFS(n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q)
				# print(dip_elem)
				dip_elem_square = (dip_elem * a0 * e) ** 2
				sum_dip_square += dip_elem_square
			# print('New f...')


	light_shift_transition = sum_dip_square * get_E_square(P, w0) / (4 * h * det_D32)
	# In MHz
	light_shift_transition /= h * 1e6
	light_shift += light_shift_transition
	print('Light Shift due to 5D3/2 in MHz: \n', light_shift_transition)
	sum_dip_square = 0

	
	'''
	Sum for all contributions from 5P3/2 ---> 7S1/2 ; 740.82 nm
	'''

	lambda_trans_7S12 = 740.82e-9
	lambda_laser = 850e-9
	nu_laser = c / lambda_laser
	nu_transition_7S12 = c / lambda_trans_7S12
	det_7S12 = nu_laser - nu_transition_7S12


	n2 = 7
	l2 = 0
	j2 = 1/2

	f_min = int(atom.I - j2)
	f_max = int(atom.I + j2)

	print('\nCalculating light shift for 5P3/2 ---> 7S1/2 contribution...')
	for f2 in range(f_min, f_max + 1, 1):
		for mf2 in range(-f2, f2 + 1, 1):
			# print('f1 = ', f1)
			if mf2 == mf1 + q and abs(f2 - f1) < 2:
				# print('mf1 = ', mf1)
				dip_elem = atom.getDipoleMatrixElementHFS(n1,l1,j1,f1,mf1,n2,l2,j2,f2,mf2,q)
				# print(dip_elem)
				dip_elem_square = (dip_elem * a0 * e) ** 2
				sum_dip_square += dip_elem_square
			# print('New f...')


	light_shift_transition = sum_dip_square * get_E_square(P, w0) / (4 * h * det_D32)
	# In MHz
	light_shift_transition /= h * 1e6
	light_shift += light_shift_transition
	print('Light Shift due to 7S1/2 in MHz: \n', light_shift_transition)
	sum_dip_square = 0

	return light_shift


def detuning(lambda_laser=850e-9):
	lambda_trans_32 = 780e-9
	lambda_trans_12 = 795e-9
	nu_laser = c / lambda_laser
	nu_transition_12 = c / lambda_trans_12
	nu_transition_32 = c / lambda_trans_32

	det_D32 = nu_laser - nu_transition_32
	det_D12 = nu_laser - nu_transition_12

	return det_D12, det_D32


def get_E_square(P, w0):
	E2 = 2 * get_I(P, w0) / (epsilon_0 * c)
	return E2


def get_I(P, w0):
	I = 2 * P / (np.pi*w0**2)
	return I


def get_ls_mfs_5S(atom, f, q, P, w0):
	mf_matrix = np.array([])
	ls_matrix = np.array([])

	for mf in range(-f, f + 1, 1):
		detuning_s_state = s_states(atom,f,mf,q, P, w0)
		mf_matrix = np.append(mf_matrix, mf)
		ls_matrix = np.append(ls_matrix, detuning_s_state)

	label = r'5S$_{1/2}$, f = ' + str(f) + r', $\sigma_{pol} = $' + str(q)
	plt.plot(mf_matrix, ls_matrix, '_', markersize=30, label=label)
	plt.ylabel('Light shift (MHz)', fontsize=fontsize)
	plt.xlabel('mf', fontsize=fontsize)
	plt.legend(fontsize=fontsize)

	return mf_matrix, ls_matrix


def get_ls_mfs_5P3_2(atom, f, q, P, w0):
	mf_matrix = np.array([])
	ls_matrix = np.array([])
	title = r'$\omega_{0} = $' + str(w0*1e6) + r'$\mu$' +'m' + ', P = ' + str(P*1e3) + ' mW'

	for mf in range(-f, f + 1, 1):
		detuning_p_32_state = p_3_2_states(atom,f,mf,q, P, w0)
		mf_matrix = np.append(mf_matrix, mf)
		ls_matrix = np.append(ls_matrix, detuning_p_32_state)

	label = r'5P$_{3/2}$, f = ' + str(f) + r', $\sigma_{pol} = $' + str(q)
	plt.plot(mf_matrix, ls_matrix, '_', markersize=30, label=label)
	plt.ylabel('Light shift (MHz)', fontsize=fontsize)
	plt.xlabel(r' Zeeman sub-level, $m_{f}$', fontsize=fontsize)
	plt.xticks(ticks=None)
	plt.legend(fontsize=fontsize)
	plt.title(title, fontsize=fontsize)

	return mf_matrix, ls_matrix


def get_ls_transition5S_5P32(atom, f_ground, mf_ground, f_excited, mf_excited, q, P, w0):
	detuning_s_state = s_states(atom, f_ground, mf_ground, q, P, w0)
	detuning_p_state = p_3_2_states(atom, f_excited, mf_excited, q, P, w0)

	shift = detuning_p_state - detuning_s_state

	return shift


def get_plot_transitions_OP_pushout(atom, w0):
	# Pushout transition
	q = 0 # trap lasers


	p = np.linspace(0, 5, 5)
	p *= 1e-3


	ls_shift_pushout = np.array([])
	ls_shift_OP = np.array([])
	plt.figure('Power light shift')
	# Pushout transition
	for power in p:
		# Pushout transition
		f_ground = 2
		mf_ground = 2
		f_excited = 3
		mf_excited = mf_ground + 1
		shift = get_ls_transition5S_5P32(atom, f_ground, mf_ground, f_excited, mf_excited, q, power, w0)
		ls_shift_pushout = np.append(ls_shift_pushout, shift)

		# OP transition
		f_excited = 2
		mf_ground = 0
		mf_excited = mf_ground + 1
		shift = get_ls_transition5S_5P32(atom, f_ground, mf_ground, f_excited, mf_excited, q, power, w0)
		ls_shift_OP = np.append(ls_shift_OP, shift)

	plt.plot(p*1e3, ls_shift_pushout, label='Pushout transition')
	plt.plot(p*1e3, ls_shift_OP, label='Optical pumping transition')
	plt.xlabel('Laser power (mW)', fontsize=fontsize)
	plt.ylabel('Light shift (MHz)', fontsize=fontsize)
	plt.xlim((0, 5))
	plt.ylim((0, 40))
	plt.legend(fontsize=fontsize)
	plt.grid(linewidth=.2)
	plt.show()


def get_ls_OP_pushout(laser_power, q=0):
	atom = Rubidium87()
	w0 = 1.1e-6

	# Quantum numbers
	f_ground = 2
	mf_ground = 2
	# Pushout transition
	f_excited = 3
	mf_excited = mf_ground + 1
	shift_pushout = get_ls_transition5S_5P32(atom, f_ground, mf_ground, f_excited, mf_excited, q, laser_power, w0)


	# OP transition
	f_excited = 2
	mf_ground = 0
	mf_excited = mf_ground + 1
	shift_OP = get_ls_transition5S_5P32(atom, f_ground, mf_ground, f_excited, mf_excited, q, laser_power, w0)


	return shift_OP, shift_pushout



