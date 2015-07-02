#!/usr/bin/env python

# This Script looks for couplings higher than a certain threshold in all the 
# coupling calculationis contained in the PWD. The PWD should contain folders
# named V_chrom1.chrom2, containing the V_chrom1.chrom2.log logfiles from
# coupling calculations with G09. At the moment, it works only with gdvh23dev
# version of G09. I will extend it to the gdvh36 version in case I will use it.
# One can select the threshold, the unit and the type of coupling to be checked
# through command line options.

# Import standard Python modules
import os
import re
import argparse as arg


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''
    Analyzes the couplings obtained from G09 EET calculations.''')

    # Optional arguments
    # Threshold
    parser.add_argument('-t', '--thresh', default=0.65, type=float, help='''
    Threshold under which couplings will be ignored.''')

    # Unit
    parser.add_argument('-u', '--unit', default='eV', choices=['eV', 'cm'],
    help='''Unit of the Threshold.''')

    # Coupling type
    parser.add_argument('-c', '--coup', default='coulomb',
    choices=['coulomb', 'tot'], help='''Unit of the Threshold.''')

    args = parser.parse_args()

    return args


# Define a function to check whether an iterable is empty or not
def empty(iterable):
    try:
        return all(map(empty, iterable))
    except TypeError:
        return False


# =========================
# The Program Starts Here
# =========================

if __name__ == '__main__':

    args = options()
    thresh = args.thresh
    unit = args.unit
    coup_type = args.coup

    # Set up the column containing the coupling in the desired unit
    if unit == 'eV':
        col = 8
    if unit == 'cm':
        col = 10

    # Set up the string to be matched for the desired coupling
    if coup_type == 'coulomb':
        search_string = 'Coulomb term'
    if coup_type == 'tot':
        search_string = 'TOTAL COUPLING'

    # Couplings defined as the folders whose name equals
    # 'V_' followed by some digits
    couplings = sorted(
        filter(lambda x: re.match('V_\d+', x), os.listdir(os.getcwd())))

    # Initialize a dictionary to store couplings higher than the threshold.
    # The key will be the name of the coupling, the value will be a list
    # containing the couplings higher than the threshold.
    suspect_data = {}

    for coupling in couplings:
        suspect_data[coupling] = []
        with open(os.path.join(coupling, '%s.log' % coupling), 'r') as f:
            for line in f:
                if search_string in line:
                    coup = float(line.split()[col])
                    if abs(coup) > thresh:
                        suspect_data[coupling].append(coup)

    if empty(suspect_data.values()):
        print('''
        No couplings higher than %f %s have been found.''' % (thresh, unit))
    else:
        print('''
        Couplings higher than %f %s have been saved
        in coup_analysis.dat''' % (thresh, unit))

        with open('coup_analysis.dat', 'w') as data:
            for k, v in suspect_data.iteritems():
                if v:
                    data.write('%s\n\n' % k)
                    for c in v:
                        data.write('\t%10.6f\n' % c)
                    data.write('\n')
