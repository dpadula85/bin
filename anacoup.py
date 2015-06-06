#!/usr/bin/env python

import sys, os, re

if len(sys.argv) > 1:
  threshold = float(sys.argv[1])
else:
  threshold = 0.65

def empty(seq):
  try:
    return all(map(empty, seq))
  except TypeError:
    return False

couplings = sorted(filter(lambda x: re.match('V_\d+', x), os.listdir(os.getcwd())))

suspect_data = {}

for coupling in couplings:
  suspect_data[coupling] = []  
  with open(os.path.join(coupling, '%s.log' % coupling), 'r') as f:
    for line in f:
      if 'Coulomb term' in line:
        coup = float(line.split()[8])
        if abs(coup) > threshold:
          suspect_data[coupling].append(coup)  

if empty(suspect_data.values()):
  print('No couplings higher than %6.4f eV have been found.' % threshold)
else:
  print('Couplings higher than %6.4f eV have been printed in strange_coup.dat' % threshold)
  with open('strange_coup.dat', 'w') as data:
    for k, v in suspect_data.iteritems():
      if v:
        data.write('%s\n\n' % k)
        for c in v:
          data.write('\t%10.6f\n' % c)
        data.write('\n')
