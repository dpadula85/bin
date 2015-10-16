#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

diag = np.loadtxt('diag.dat')

# Get coefficients from the diagonalized matrix
coeff = diag[:diag.shape[0]/2,2:]

fig, ax = plt.subplots()
im = ax.pcolor(np.square(coeff))
fig.colorbar(im)

# Generate labels
xlabels = np.arange(1, len(coeff[0])+1)

# X axis
# Hide main labels and assign to minor labels their value
plt.xlabel("Local States")
ax.set_xticks(xlabels)
ax.set_xticklabels(xlabels, visible=False)
ax.set_xticks(xlabels - 0.5, minor=True)
ax.set_xticklabels(xlabels, minor=True)
ax.set_xlim(right=xlabels[-1])


# Y axis
# Hide main labels and assign to minor labels their value
plt.ylabel("Excitonic States")
ax.set_yticks(xlabels)
ax.set_yticklabels(xlabels, visible=False)
ax.set_yticks(xlabels - 0.5, minor=True)
ax.set_yticklabels(xlabels, minor=True)
ax.set_ylim(top=xlabels[-1])

plt.setp(ax.xaxis.get_minorticklabels(), rotation=-45 )
plt.setp(ax.yaxis.get_minorticklabels(), rotation=-45 )

plt.grid()
plt.show()
