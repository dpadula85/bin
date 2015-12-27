#!/usr/bin/env python

import sys
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

    parser.add_argument('--fit', choices=['single', 'global'], default='single',
    help='''Type of fit to use. Single : each peak fitted alone, Global : whole data set fitted; the components
    are also plotted.''')

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


def gaussians(x, *parms):

    y = np.zeros_like(x)

    for i in range(0, len(parms), 3):
        A = parms[i]
        avg = parms[i + 1]
        sigma = parms[i + 2]
        y += A * np.exp( -(x - avg)**2/ (2. * sigma**2))

    return y


def lorentzians(x, *parms):

    y = np.zeros_like(x)

    for i in range(0, len(parms), 3):
        A = parms[i]
        avg = parms[i + 1]
        gamma = parms[i + 2]
        y += A * (1 / np.pi) * gamma / (gamma**2 + (x - avg)**2)

    return y


if __name__ == '__main__':

    args = options()

    # acc = args.accuracy
    lineshape = args.lineshape
    data = np.loadtxt(args.filename)

    # Get DataSet and add it to the plot
    x = data[:,0]
    y = data[:,1]
    title = "%s fit" % args.fit
    plt.title(title.title())
    plt.plot(x, y, lw=2, color='black', label='Data Set')

    # Find maxima and minima in the dataset
    peaks = findpeaks(x, y)

    print
    print(banner("Deconvolution", "=", 60))
    print
    print(" > Found %d peaks:" % len(peaks))
    print

    totguess = np.array([])
    for i, p in enumerate(peaks, start=1):

        xm = p[0]
        ym = p[1]

        # Set the initial guess for the fitting procedure
        A = 50.0
        avg = xm
        wid = 3
        parms = np.array([A, avg, wid])
        totguess = np.r_[totguess, parms]

        # Add points for the maxima and minima to the plot
        marker = plt.plot(xm, ym, 'o', ms=8, label='Det. Max. %d' % i)
        col = marker[0].get_color()

        if args.fit == 'single':

            if lineshape == 'gau':
                funct = gaussians

            elif lineshape == 'lor':
                funct = lorentzians

            # Fit the single maximum with a function
            popt, pcov = curve_fit(funct, x, y, p0=parms)
            y_fitted = funct(x, *popt)

            # Summary of the fitting procedure
            print(banner("Peak %d" % i, "=", 30))
            print(" > Area  : %10.2f" % np.abs(popt[0]))
            print(" > Max   : %10.2f" % popt[1])
            print(" > Sigma : %10.2f" % popt[2])
            print

            # Plot the result of the fitting
            peak = plt.plot(x, y_fitted, color=col, lw=2, linestyle='dashed', label='Deconv. Peak %d' % i)

    # Here we're out of the for cycle! Do operations for the global fit
    if args.fit == 'global':

        # Fit the data set with a sum of n functions, where n is the
        # number of peaks found by the findpeaks function
        if lineshape == 'gau':

            try:
                popt, pcov = curve_fit(gaussians, x, y, p0=totguess)

            except RuntimeError:
                print(banner("ERROR", "=", 60))
                print(" Curve fitting failed!")
                sys.exit()

            y_totfit = gaussians(x, *popt)

        elif lineshape == 'lor':

            try:
                popt, pcov = curve_fit(lorentzians, x, y, p0=totguess)

            except RuntimeError:
                print(banner("ERROR", "=", 60))
                print(" Curve fitting failed!")
                sys.exit()

            y_totfit = lorentzians(x, *popt)

        # Plot the function that fits experimental data
        tot = plt.plot(x, y_totfit, color='black', lw=3, linestyle='dashed', label='Total fit')

        plt.gca().set_color_cycle(None)

        if lineshape == 'gau':
            funct = gaussians

        elif lineshape == 'lor':
            funct = lorentzians

        # Plot the single functions constituting the sum function that fits
        # the experimental data set
        j = 1
        for i in range(0, len(totguess), 3):

            A = popt[i]
            avg = popt[i + 1]
            wid = popt[i + 2]
            parms = np.array([A, avg, wid])

            y_fitted = funct(x, *parms)
            peak = plt.plot(x, y_fitted, lw=2, linestyle='dashed', label='Component %d' % j)

            # Summary of the fitting procedure
            print(banner("Component %d" % j, "=", 30))
            print(" > Area  : %10.2f" % np.abs(parms[0]))
            print(" > Max   : %10.2f" % parms[1])
            print(" > Sigma : %10.2f" % parms[2])
            print

            j += 1

    if len(peaks) > 0:
        plt.legend().draw_frame(False)
        plt.show()
