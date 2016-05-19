#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

ev2wn = 8065.73

# Define Harmonic Diabatic States
def H_GS_m1(w, q):
    return 0.5 * w * q**2

def H_ES_m1(E0, w, q0, q):
    return E0 + 0.5 * w * (q - q0)**2

def H_GS_diab(w, q1, q2):
    return H_GS_m1(w, q1) + H_GS_m1(w, q2)

def H_ES_diab(E0, w, q0, q1, q2):
    return H_GS_m1(w, q1) + H_ES_m1(E0, w, q0, q2)

def H_ES_ad1(E0, w, q0, q1, q2):
    return (H_ES_diab(E0, w, q0, q1, q2) + H_ES_diab(E0, w, q0, q2, q1)) / 2 - np.sqrt(0.25 * (H_ES_diab(E0, w, q0, q1, q2) - H_ES_diab(E0, w, q0, q2, q1))**2 + V**2)

def H_ES_ad2(E0, w, q0, q1, q2):
    return (H_ES_diab(E0, w, q0, q1, q2) + H_ES_diab(E0, w, q0, q2, q1)) / 2 + np.sqrt(0.25 * (H_ES_diab(E0, w, q0, q1, q2) - H_ES_diab(E0, w, q0, q2, q1))**2 + V**2)
#
# Define parameters in wavenumbers
#
linterm = 2244 #/ ev2wn
V = -340 #/ ev2wn
w = 1501 #/ev2wn
x0 = linterm / np.sqrt(w)**3
delta0 = np.sqrt(w) * x0
E0 = 29600 #/ ev2wn
Er = 0.5 * w * delta0**2

q1 = q2 = np.linspace(-2*delta0, 2*delta0, 50)
x, y = np.meshgrid(q1, q2)

#
# Options
#
fig = plt.figure(figsize=(30,30))
gs = gridspec.GridSpec(2, 2)
gs.update(wspace=0.05, hspace=0.05)

#
# Plot Diabatic states PESs
#
ax = plt.subplot(gs[0, 0],projection='3d')
ax.set_title('Diabatic States')
ax.set_xlabel(r'q$_1$')
ax.set_ylabel(r'q$_2$')
ax.set_zlabel(r'E (eV)')

# GS
z = H_GS_diab(w, x, y) + 0.5*E0
z = z.reshape(x.shape)
# ax.plot_surface(x, y, z, color='red', alpha=0.5, label='GS')

# L1 state
z1 = H_ES_diab(E0, w, delta0, x, y)
z1 = z1.reshape(x.shape)
ax.plot_surface(x, y, z1, color='lime', alpha=0.5, label='GS')

# L2 state
z2 = H_ES_diab(E0, w, delta0, y, x)
z2 = z2.reshape(x.shape)
ax.plot_surface(x, y, z2, color='blue', alpha=0.5, label='GS')

#
# Plot Adiabatic states PESs
#
ax1 = plt.subplot(gs[0, 1],projection='3d')
ax1.set_title('Adiabatic States')
ax1.set_xlabel(r'q$_1$')
ax1.set_ylabel(r'q$_2$')
ax1.set_zlabel(r'E (eV)')

# GS
z = H_GS_diab(w, x, y) + 0.5*E0
z = z.reshape(x.shape)
# ax1.plot_surface(x, y, z, color='red', alpha=0.5, label='GS')

# A1 state
z3 = H_ES_ad1(E0, w, delta0, x, y)
z3 = z3.reshape(x.shape)
ax1.plot_surface(x, y, z3, color='cyan', alpha=0.5, label='GS')

# A2 state
z4 = H_ES_ad2(E0, w, delta0, x, y)
z4 = z4.reshape(x.shape)
ax1.plot_surface(x, y, z4, color='orange', alpha=0.5, label='GS')

#
# Plot all surfaces along q- coordinate
#
ax2 = plt.subplot(gs[1, :])
ax2.set_xlabel(r'q$_-$')
ax2.set_ylabel(r'E (cm$^{-1}$)')
ax2.set_ylim(E0, E0 + 0.25*E0)

q_min = (q1 + q2) / np.sqrt(2)

# L1 state
z5 = H_ES_diab(E0, w, delta0, -q1, q1)
ax2.plot(q_min, z5, color='lime', lw=3, ls='dashed', label='L1')

# L2 state
z6 = H_ES_diab(E0, w, delta0, q1, -q1)
ax2.plot(q_min, z6, color='blue', lw=3, ls='dashed', label='L2')

# A1 state
z7 = H_ES_ad1(E0, w, delta0, q1, -q1)
ax2.plot(q_min, z7, color='cyan', lw=2, label='A1')

# A2 state
z8 = H_ES_ad2(E0, w, delta0, q1, -q1)
ax2.plot(q_min, z8, color='orange', lw=2, label='A2')

# plt.savefig('diab.svg', dpi=1200, transparent=True, size=(30,30))
plt.show()
