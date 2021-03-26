import numpy as np
from arc import *
from scipy.constants import h, hbar, k, c, epsilon_0, e
import matplotlib.pyplot as plt


'''
Some plot settings
'''
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

# TeX settings
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}']  # for \text command



data_cal_v_power = np.loadtxt('VtoP_data.dat')

plt.plot(data_cal_v_power[:, 0], data_cal_v_power[:, 1], 'o', markerfacecolor='blue')

plt.show()

data_cal_v_ls = np.loadtxt('VtoLS_data.dat')

plt.plot(data_cal_v_ls[:, 0], data_cal_v_ls[:, 1], 'o', markerfacecolor='blue')
plt.show()