#!/usr/bin/env python

import numpy as np
import argparse as arg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import signal


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Deconvolutes Experimental Data in separated Peaks.')

    # Optional arguments
    parser.add_argument('-f', '--filename', help='''File data.dat from mdanalyzer.''')

    # parser.add_argument('-a', '--accuracy', type=float, default=10., help='''Tolerance
    # interval for the detection of another maximum/minimum along the x axis.''')

    parser.add_argument('-ls', '--lineshape', choices=['gau', 'lor'], default='gau',
    help='''Type of function for the deconvolution of experimental peaks.''')

    args = parser.parse_args()

    return args


def banner(text=None, ch='=', length=78):
    """Return a banner line centering the given text.
    
        "text" is the text to show in the banner. None can be given to have
            no text.
        "ch" (optional, default '=') is the banner line character (can
            also be a short string to repeat).
        "length" (optional, default 78) is the length of banner to make.

    Examples:
        >>> banner("Peggy Sue")
        '================================= Peggy Sue =================================='
        >>> banner("Peggy Sue", ch='-', length=50)
        '------------------- Peggy Sue --------------------'
        >>> banner("Pretty pretty pretty pretty Peggy Sue", length=40)
        'Pretty pretty pretty pretty Peggy Sue'
    """
    if text is None:
        return ch * length

    elif len(text) + 2 + len(ch)*2 > length:
        # Not enough space for even one line char (plus space) around text.
        return text

    else:
        remain = length - (len(text) + 2)
        prefix_len = remain / 2
        suffix_len = remain - prefix_len
    
        if len(ch) == 1:
            prefix = ch * prefix_len
            suffix = ch * suffix_len

        else:
            prefix = ch * (prefix_len/len(ch)) + ch[:prefix_len%len(ch)]
            suffix = ch * (suffix_len/len(ch)) + ch[:suffix_len%len(ch)]

        return prefix + ' ' + text + ' ' + suffix


def findpeaks(x, y):

    _max = signal.argrelmax(y, order=5)[0]
    _min = signal.argrelmin(y, order=5)[0]
    maxs = [[x[p], y[p]] for p in _max]
    mins = [[x[p], y[p]] for p in _min]

    peaks = maxs + mins
    peaks.sort(key=lambda x: x[0])

    return peaks


def gaussian_peak(x, *parms):

    A, avg, sigma = parms

    return A * np.exp( -(x - avg)**2 / (2. * sigma**2))


def lorentzian_peak(x, *parms):

    A, avg, gamma = parms

    return A * (1 / np.pi) * gamma / (gamma**2 + (x - avg)**2) 


def fit_peak(x, y, guess, function='gau'):
    
    if function == 'gau':
        function = gaussian_peak

    elif function == 'lor':
        function = lorentzian_peak

    popt, pcov = curve_fit(function, x, y, p0=guess)
    y_fitted = function(x, *popt)

    return popt, y_fitted


if __name__ == '__main__':

    args = options()

    # acc = args.accuracy
    lineshape = args.lineshape
    data = np.loadtxt(args.filename)

    # Get DataSet and add it to the plot
    x = data[:,0]
    y = data[:,1]
    plt.plot(x, y, lw=2, color='black', label='Data Set')

    # Find maxima and minima in the dataset
    peaks = findpeaks(x, y)

    print
    print(banner("Deconvolution", "=", 60))
    print
    print(" > Found %d peaks:" % len(peaks))
    print

    for i, p in enumerate(peaks, start=1):

        xm = p[0]
        ym = p[1]

        # Fit the peak with a function
        A = 50.0
        avg = xm
        sigma = 3
        parms = np.array([A, avg, sigma])
        popt, y_fitted = fit_peak(x, y, parms, lineshape)

        # Summary of the fitting procedure
        print(banner("Peak %d" % i, "=", 30))
        print(" > Area  : %10.2f" % np.abs(popt[0]))
        print(" > Max   : %10.2f" % popt[1])
        print(" > Sigma : %10.2f" % popt[2])
        print

        # Add data to the plot
        marker = plt.plot(xm, ym, 'o', ms=8, label='Det. Max. %d' % i)
        col = marker[0].get_color()
        peak = plt.plot(x, y_fitted, color=col, lw=2, linestyle='dashed', label='Deconv. Peak %d' % i)


    if len(peaks) > 0:
        plt.legend().draw_frame(False)
        plt.show()
    pass
