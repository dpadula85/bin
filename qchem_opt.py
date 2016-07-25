#!/usr/bin/env python

import sys

def skiplines(openfile, nlines=0):
    '''Skips nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.'''

    for i in range(nlines):
        next(openfile)

    return next(openfile)

fname = sys.argv[1]
opt = []
struct = []

with open(fname) as f:
    for line in f:
        if "Optimization Cycle" in line:
            line = skiplines(f, 3)
            data = line.split()
    
            while len(data) == 5:
                data[2:] = map(float, data[2:])
                struct.append(data[1:])
                line = skiplines(f)
                data = line.split()

            opt.append(struct)
            struct = []

for i, step in enumerate(opt, start=1):
    print len(step)
    print "Optimization Step %d" % i

    for atom in step:
        print "%3s %14.8f %14.8f %14.8f" % tuple(atom)
