import numpy as np
import matplotlib.pyplot as plt
from arc import *
import ls_functions as ls
import sys

fontsize = 16

atom = Rubidium87()

q = int(0) # Photon Polarization
laser_power = 5e-3
P = laser_power
w0 = 1.1e-6
B = 6.5
ls.plot_pushout_w_ls_zeeman(atom, B, laser_power, w0=1.1e-6)
ls.plot_op_w_ls_zeeman(atom, B, laser_power, w0=1.1e-6)

sys.exit(0)

'''

Here I calculate the light shift with respect to the HFS
on the mf levels of a particular f state. Will plot shift againt mf

'''
f = 2
mf_5S12_matrix, ls_5S12_matrix = ls.get_ls_mfs_5S(atom, f, q, P, w0)
f = 3
mf_5P32_matrix, ls_5P32_matrix = ls.get_ls_mfs_5P3_2(atom, f, q, P, w0)
f = 2
mf_5P32_matrix, ls_5P32_matrix = ls.get_ls_mfs_5P3_2(atom, f, q, P, w0)

# q = 1 # Photon Polarization
# f = 2
# mf_5S12_matrix, ls_5S12_matrix = ls.get_ls_mfs_5S(atom, f, q, P, w0)
# f = 2
# mf_5P32_matrix, ls_5P32_matrix = ls.get_ls_mfs_5P3_2(atom, f, q, P, w0)

plt.grid(linewidth=.2)



'''

Below I calculate the effective shift on the relevant transitions

Pushout: 5S1/2; F, mf = 2,2 ---> 5P3/2; F, mf = 3, 2, assuming q=0 (linear pol trpa)
optical pumping: 5S1/2; F, mf = 2, vaious ---> 5P3/2; F, mf = 2, various + 1, assuming q=0 (linear pol trap)
'''
ls.get_plot_transitions_OP_pushout(atom, w0)

ls_OP, ls_pushout = ls.get_ls_OP_pushout(4e-3)
print(ls_pushout)
print(ls_OP)


