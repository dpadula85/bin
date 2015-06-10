#!/usr/bin/env python

#Extracts data from results.out, spec.OD.dat and spec.CD.dat and prepares
#a table for faster plotting with qtiplot.

energies = []
f_osc = []
rot = []

with open('results.out', 'r') as f:
    for line in f:
        energies.append(float(line.split()[1]))
        f_osc.append(float(line.split()[2]))
        rot.append(float(line.split()[-1]))
    
        
wavelengths = []
UV = []
CD = []

with open('spec.OD.dat', 'r') as f:
    for line in f:
        if line.startswith('#'):
            pass
        else:
            wavelengths.append(float(line.split()[0]))
            UV.append(float(line.split()[1]))
            
with open('spec.CD.dat', 'r') as f:
    for line in f:
        if line.startswith('#'):
            pass
        else:
            CD.append(float(line.split()[1]))
            
f = open('results_qti.txt', 'w')

for i in range(len(energies)):
    f.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (wavelengths[i], UV[i], CD[i], 1240 / energies[i], f_osc[i], f_osc[i] * 50000, rot[i], rot[i] / 3))
    
for j in range(i + 1, len(wavelengths)):
    f.write('%f\t%f\t%f\n' % (wavelengths[j], UV[j], CD[j]))
    
f.close()
