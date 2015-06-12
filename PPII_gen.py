#!/usr/bin/env python

# Generates a series of files to be used in PyMol to build polypeptides in
# polyproline II conformation

# import standard Python Modules
from numpy import arange

# Definition of the zone of Ramachadran plot describing PPII
# Angles in degrees
phi_min = -80.0
phi_max = -60.0

psi_min = +140.0
psi_max = +180.0

# Step in degrees
step = 5.0

# Aminoacid
AA = 'A'

# Number of residues
nres = 21

# List of residues
# Potentially it can be extended to the case where each residue is different
res = [AA] * nres

for phi in arange(phi_min, phi_max + 1, step):
    for psi in arange(psi_min, psi_max + 1, step):
        with open('seq_%s_%s.txt' %
                  (int(abs(phi)), int(abs(psi))), 'w') as seq_file:
            for i in range(len(res)):
                seq_file.write('%s\t%5.1f\t%5.1f\n' % (res[i], phi, psi))
            seq_file.close()
