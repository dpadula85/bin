#!/usr/bin/env python

import os
import numpy as np
import argparse as arg
import matplotlib.pyplot as plt

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
    parser.add_argument('-f', '--filename', default='data.dat', help='''File data.dat from mdanalyzer.''')

    parser.add_argument('--c1', default=1, type=int, help='''Time Column.''')

    parser.add_argument('--c2', default=['2'], nargs='+', help='''Data Columns.''')

    parser.add_argument('-tu', '--timeunit', choices=["s", "ms", "mus", "ns", "ps", "fs"],
                        type=str, default="ps", help='''Unit of the time series''')

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
    # Time step in seconds, assuming points are evenly spaced
    #
    ts = (x[1] - x[0]) * factor

    #
    # Calculate freqs in Hz and convert to wavenumbers
    #
    xf = np.fft.fftfreq(N, d=ts) / c

    #
    # Calculate the FFT of y and normalise
    #
    yf = (2.0 / N) * np.fft.fft(y)

    #
    # Return only the positive half of the freqs and the real part of the FFT
    # i.e. the Cosine Transform
    #
    xf = xf[:N/2]
    yf = np.abs(yf[:N/2].real)

    return xf, yf


if __name__ == '__main__':


    args = options()
    f = args.filename
    c1 = args.c1 - 1
    c2 = map(lambda x: x - 1, extend_compact_list(args.c2))
    time_factor = time_units[args.timeunit]

    data = np.loadtxt(f)
    x = data[:,c1]
    y = data[:,c2]
    acf_y = acf(y)

    #
    # Get the Cosine Transform of the autocorrelation function
    #
    freqs, specden_ft_part = cos_transf(x, acf_y, factor=time_factor)

    #
    # Calculate the total Spectral Density and convert to wavenumbers
    #
    prefac = freqs / (k_B * T * np.pi)
    specden = prefac * specden_ft_part * eV2wn

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

    plt.plot(freqs, specden)
    plt.show()
