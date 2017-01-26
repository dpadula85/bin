#!/usr/bin/env python

import sys
import numpy as np
import argparse as arg
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Smoothes data in the input file applying a Savitzky-Golay filter.')

    # Optional arguments
    parser.add_argument('-f', '--filename', help='''Input File.''')

    parser.add_argument('--show', help='''Show the plot in an external window.''',
    default=False, action='store_true')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = options()

    data = np.loadtxt(args.filename)
    basename = '.'.join(args.filename.split('.')[:-1])

    # Get DataSet and add it to the plot
    x = data[:,0]
    y = data[:,1]
    yhat = savgol_filter(y, 51, 3)
    hdr = "\n Smoothed Data\n"
    np.savetxt("%s.smooth.out" % basename, np.c_[x, yhat], fmt="%18.6e", header=hdr)

    plt.plot(x, y, color="b", label="Data Set")
    plt.plot(x, yhat, color="k", label="Smoothed")

    if args.show:
        plt.show()
