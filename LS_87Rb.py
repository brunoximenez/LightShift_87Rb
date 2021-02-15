import numpy as np
import matplotlib.pyplot as plt
from arc import *
import light_shift_calc as ls
import sys

color_edge_list = [(0.,0.4,0.4),(1,0.2,0.2),(0.5,0,0.5),(0.2,0.2,1),(0.2,0.2,0.2),(1,0.6,0.),(0,0.6,0),(0.2,0.2,0),(0.4,0.4,1)]
color_face_list = [(0.,0.8,0.8),(1,0.7,0.7),(0.9,0.2,0.9),(0.7,0.7,1),(0.7,0.7,0.7),(1,0.8,0.4),(0,0.9,0),(0.5,0.5,0),(0.6,0.6,1)]

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

element = Rubidium87()

laser_power = np.linspace(1e-3, 5e-3, 5)
w0 = 1.1e-6


ls_d1 = np.array([])
ls_d2 = np.array([])
u0 = np.array([])


det_12, det_32 = ls.detuning(850e-9)


for p in laser_power:
    light_shift_d1, light_shift_d2, trap_depth = ls.light_shift_Ds(element, p, w0, det_12, det_32)
    ls_d1 = np.append(ls_d1, light_shift_d1)
    ls_d2 = np.append(ls_d2, light_shift_d2)
    u0 = np.append(u0, trap_depth)


title = 'w0 = ' + str(w0)
plt.subplot(211)
plt.title(title)
plt.plot(u0, ls_d1, label='D1 transition')
plt.plot(u0, ls_d2, label='D2 transition')
plt.ylabel('Light shift (MHz)')
plt.xlim((np.amin(u0), np.amax(u0)))
plt.legend()
plt.grid()
plt.subplot(212)
plt.plot(u0, laser_power*1e3)
plt.ylabel('Laser power (mW)')
plt.xlabel('Trap depth (MHz)')
plt.xlim((np.amin(u0), np.amax(u0)))
plt.ylim((np.amin(laser_power*1e3), np.amax(laser_power*1e3)))
plt.grid()
plt.show()
sys.exit(0)

#Find the different lightshift made by the 1005 on the 6P_3/2, m_J = 3/2

n_list = np.arange(4,80,1)


fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(111)

i  = 0
for l in [0,2]:
    if l == 0 :
        j = 0.5 
        label = r'$6P_{1/2} \rightarrow nS_{1/2}$'

    else : 
        j = 1.5
        label = r'$6P_{1/2} \rightarrow nD_{3/2}$'
        
    
    #for different polarization
    for polar in [-1,0,1]:
        
        rabi_list = []
        det_list = []
        lightshift_list = []
    
        for n in n_list :
            state = [n,l,j]
            
            rabi = atom.getRabiFrequency(
                n1 = 6,
                l1 = 1,
                j1 = 0.5,
                mj1 = 0.5,
                n2 = state[0],
                l2 = state[1],
                j2 = state[2],
                q = polar,
                laserPower = 150e-3,
                laserWaist = 3.4e-6
                )
            
            transition_freq = atom.getTransitionFrequency(
                n1 = 6,
                l1 = 1,
                j1 = 0.5,
                n2 = state[0],
                l2 = state[1],
                j2 = state[2]
                )
        
            detuning_THz = (298.145*1e12-np.abs(transition_freq))
            lightshift = rabi**2/4./(detuning_THz*2*np.pi)
            lightshift_MHz = lightshift/2/np.pi*1e-6
            
            rabi_list.append(rabi)
            det_list.append(detuning_THz)
            lightshift_list.append(lightshift_MHz)
        
        lightshift_list = np.asarray(lightshift_list)
        sum_lightshift = np.sum(lightshift_list[:55])
        
        ax.plot(
            n_list, 
            lightshift_list, 
            marker = 'o',
            markeredgewidth = 2,
            markersize = 8,
            markerfacecolor = color_face_list[i],
            markeredgecolor = color_edge_list[i],
            linestyle = '',
            label = label + ' polar = %d'%polar + ' total = %.1f MHz'%sum_lightshift
            )
        i += 1

ax.legend(loc = 1, numpoints = 1, fontsize = 14).get_frame().set_alpha(0.8)
ax.set_xlabel('principal quantum number n')
ax.set_ylabel('Lightshift (MHz)')
#plt.savefig('lightshift_1004')
plt.show()
#
#print("6 P_{3/2} -> %d D_{%.1f}"%(state[0],state[2]))
#print('rabi : %.1f GHz'%(rabi/2/np.pi*1e-6))
#print('detuning : %.2f THz'%(detuning_THz*1e-12))
#print('lightshift (150mW) : %.1f MHz'%lightshift_MHz)