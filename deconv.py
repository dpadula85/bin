#!/usr/bin/env python

import numpy as np
import argparse as arg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import peakdetect as pd
import util as u

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Deconvolutes Experimental Data in Gaussian Peaks.')

    # Optional arguments
    parser.add_argument('-f', '--filename', help='''File data.dat from mdanalyzer.''')

    args = parser.parse_args()

    return args


def gaussian_peak(x, *parms):

    A, avg, sigma = parms

    return A * np.exp(-(x-avg)**2/(2.*sigma**2))


def fit_peak(x, y, guess):
    
    popt, pcov = curve_fit(gaussian_peak, x, y, p0=guess)
    y_fitted = gaussian_peak(x, *popt)

    return popt, y_fitted


if __name__ == '__main__':

    args = options()

    data = np.loadtxt(args.filename)
    x = data[:,0]
    y = data[:,1]

    # Find maxima and minima in the dataset
    _max, _min = pd.peakdetect(y, x, 10, 0.30)
    peaks = _max + _min
    peaks.sort(key=lambda x: x[0])

    print(u.banner("Deconvolution", "=", 60))
    print
    print(" > Found %d peaks:" % len(peaks))
    print

    result = {}
    for i, p in enumerate(peaks, start=1):

        xm = p[0]
        ym = p[1]
        A = 50.0
        avg = xm
        sigma = 3
        parms = np.array([A, avg, sigma])
        popt, y_fitted = fit_peak(x, y, parms)
        result[i] = y_fitted
        print(u.banner("Peak %d" % i, "=", 30))
        print(" > Area  : %10.2f" % popt[0])
        print(" > Max   : %10.2f" % popt[1])
        print(" > Sigma : %10.2f" % popt[2])
        print

    plt.plot(x, y, lw=2, color='black', label='Data Set')
    for i, peak in result.iteritems():

        plt.plot(x, peak, lw=2, linestyle='dashed', label='Peak %d' % i)

    plt.legend()
    plt.show()
    pass
