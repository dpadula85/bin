#!/usr/bin/env python

import os
import numpy as np
import argparse as arg
import matplotlib.pyplot as plt
from matplotlib import gridspec

eV2wn = 8065.544005

# Speed of light in cm / s
c = 2.99792458e10

# k_B : Boltzmann Constant in eV / K
k_B = 8.6173324e-05

# Temp in Kelvin
T = 300.0

time_units = {"s":  1,
              "ms": 1e-3,
              "mus": 1e-6,
              "ns": 1e-9,
              "ps": 1e-12,
              "fs": 1e-15}


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''Calculates the Fast Fourier
                                Transform of the autocorrelation function
                                of a time series.''',
                                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('-f', '--filename', default='data.dat', help='''Either a time series or an acf file.''')

    parser.add_argument('--filetype', default='series', choices=['series', 'acf'],
                        help='''Type of the input file.''')

    parser.add_argument('--c1', default=1, type=int, help='''Time Column.''')

    parser.add_argument('--c2', default=['2'], nargs='+', help='''Data Columns.''')

    parser.add_argument('-tu', '--timeunit', choices=["s", "ms", "mus", "ns", "ps", "fs"],
                        type=str, default="ps", help='''Unit of the time series''')

    parser.add_argument('-bl', '--baseline', default=None, type=float,
                        help='''Starting frequency to determine the baseline.''')

    parser.add_argument('--show', help='''Show the plot in an external window.''',
    default=False, action='store_true')

    parser.add_argument('--outtype', default='ft', choices=['specden', 'ft'],
                        help='''Type of the input file.''')

    parser.add_argument('-o', '--output', type=str, default=None, help='''Output File.''')

    args = parser.parse_args()

    return args


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    # idxs = idxs.split()

    for idx in idxs:

        to_extend = idx.split('-')

        if len(to_extend) > 1:

            sel =  map(int, to_extend)
            extended += range(sel[0],sel[1]+1,1)

        else:
        
            extended.append(int(idx))

    return extended


def acf(series):
    '''Returns the autocorrelation function of a time series.'''

    N = len(series)
    avg = np.mean(series)
    c0 = np.sum((series - avg)**2) / N

    def r(j):
        return np.sum((series[:N - j] - avg) * (series[j:] - avg)) / (N - j)

    t = np.arange(N)
    acf_t = map(r, t)

    return acf_t #/ c0


def cos_transf(x, y, factor=1):
    '''Returns the Cosine Transform of y and the frequencies'''

    #
    # Number of points
    #
    N = len(x)

    #
    # Generate 4N array for the DCT
    # Leave even elements at 0, and fill odd ones with y and y reversed
    #
    dummy_y = np.zeros(4 * N)
    dummy_y[1:2*N:2] = y
    dummy_y[2*N+1::2] = y[::-1]

    #
    # Time step in seconds, assuming points are evenly spaced
    #
    ts = (x[1] - x[0]) * factor

    #
    # Calculate freqs in Hz and convert to wavenumbers
    #
    xf = np.fft.fftfreq(4 * N, d=ts / 2.0) / c
    xf = xf[:N]
    dxf = xf[1] - xf[0]

    #
    # Calculate the FFT of y and normalise
    #
    yf = (1.0 / (N * dxf)) * np.fft.fft(dummy_y).real[:N]

    return xf, yf


if __name__ == '__main__':


    args = options()
    f = args.filename
    c1 = args.c1 - 1
    c2 = map(lambda x: x - 1, extend_compact_list(args.c2))
    time_factor = time_units[args.timeunit]
    bl = args.baseline

    data = np.loadtxt(f)

    if args.filetype == 'series':
        x = data[:,c1]
        y = data[:,c2]
        acf_y = acf(y)

    elif args.filetype == 'acf':
        x = data[:,c1]
        y = None
        acf_y = data[:,c2].flatten()

    #
    # Save to a file
    #
    data = np.c_[x, acf_y]
    if not args.output:
        basename = '.'.join(f.split('.')[:-1])
        outfile = basename + '.acf.out'

    else:
        outfile = args.output + '.acf.out'

    header = "\n      Time (%s)         ACF (eV^2)\n" % args.timeunit
    np.savetxt(outfile, data, fmt="%18.6e", header=header)

    #
    # Get the Cosine Transform of the autocorrelation function
    #
    freqs, specden_ft_part = cos_transf(x, acf_y, factor=time_factor)

    #
    # Subtract a baseline from the Cosine Transform before multiplication
    # with the prefactor (see JPCB, 2013, 117, 7157)
    #
    if bl:
        idxs = np.where(freqs > bl)
        avg = np.mean(specden_ft_part[idxs])
        specden_ft_part -= avg
        specden_ft_part[specden_ft_part < 0] = 0

    #
    # Calculate the total Spectral Density and convert to wavenumbers
    #
    if args.outtype == "specden":
        prefac = freqs / (k_B * T * np.pi)
        specden = prefac * specden_ft_part * eV2wn

    elif args.outtype == "ft":
        specden = specden_ft_part

    #
    # Save to a file
    #
    data = np.c_[freqs, specden]
    if not args.output:
        basename = '.'.join(f.split('.')[:-1])
        outfile = basename + '.specden.out'

    else:
        outfile = args.output + '.specden.out'

    header = "\n Frequency (cm^-1) Spectral Density (cm^-1)\n"
    np.savetxt(outfile, data, fmt="%18.6e", header=header)

    if args.show:

        fig = plt.figure(figsize=(11.69, 8.27))

        if y is not None:
            gs = gridspec.GridSpec(3, 1)
        else:
            gs = gridspec.GridSpec(2, 1)

        #
        # Plot FT
        #
        ax = plt.subplot(gs[0])
        ax.plot(freqs, specden, label="FT")
        ax.set_xlabel(r'$\omega$ (cm$^{-1}$)', size=22)
        ax.set_ylabel(r'$J(\omega)$ (cm$^{-1}$)', size=22)

        #
        # Plot ACF
        #
        ax0 = plt.subplot(gs[1])
        ax0.plot(x, acf_y, label="ACF")
        ax0.set_xlabel('Time (%s)' % args.timeunit, size=22)
        ax0.set_ylabel('ACF (eV$^2$)', size=22)

        if y is not None:

            #
            # Plot time series
            #
            ax1 = plt.subplot(gs[2])
            ax1.plot(x, y)
            ax1.set_xlabel('Time (%s)' % args.timeunit, size=22)
            ax1.set_ylabel('E (eV)', size=22)

        plt.tight_layout()
        plt.show()
