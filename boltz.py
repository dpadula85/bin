#!/usr/bin/env python

import sys
import numpy as np
import argparse as arg

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description="Calculates Boltzmann's Populations given the energies", formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('-f', '--filename', default='ene.dat', help='''File containing the energies in a.u.''')

    parser.add_argument('-t', '--temp', default=298.15, type=float, help='''Temperature in K.''')

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args = options()
    T = args.temp

    # Boltzmann constant in Kcal/mol*K
    k = 1.987204118e-3
    au2kcalmol = 627.509469
    data = np.loadtxt(args.filename)
    DeltaE = (data - data.min()) * au2kcalmol
    pop = np.exp(-DeltaE / (k*T)) / np.sum(np.exp(-DeltaE / (k*T)))
    results = np.c_[data, DeltaE, pop*100]
    header = "\n%11s %17s %9s\n" % ("Energy (au)", "DeltaE (Kcal/mol)", "Pop. (%)")
    np.savetxt(sys.stdout, results, header=header, fmt="%13.5f %17.4f %9.2f")
