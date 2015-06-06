#!/usr/bin/env python

import sys, os, re

couplings = sorted(filter(lambda x: re.match('V_\d+', x), os.listdir(os.getcwd())))

data = open('strange_coup.dat', 'w')

for coupling in couplings:
  with open(os.path.join(coupling, '%s.log' % coupling), 'r') as f:
    for line in f:
      if 'Coulomb term' in line:
        coup = float(line.split()[8])
        if abs(coup) > 0.65:
          data.write('%s %10.4f\n' % (coupling, coup))

data.close()
