#!/usr/bin/env python

import os
import numpy as np
import argparse as arg
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import ticker, gridspec
from matplotlib import rc
import matplotlib as mpl

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots
    "font.sans-serif": [],              # to inherit fonts from the document
    "font.monospace": [],
    "pgf.preamble": [
            r"\usepackage[utf8x]{inputenc}",
            r"\usepackage[T1]{fontenc}",
            r"\usepackage{siunitx}",
        ]
    }

mpl.rcParams.update(pgf_with_latex)



eV2wn = 8065.544005

# Speed of light in cm / s
c = 2.99792458e10

# k_B : Boltzmann Constant in eV / K
k_B = 8.6173324e-05

# Temp in Kelvin
T = 300.0

time_units = {
        "s":  1,
        "ms": 1e-3,
        "mus": 1e-6,
        "ns": 1e-9,
        "ps": 1e-12,
        "fs": 1e-15
    }


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''Calculates the Fast Fourier
                                Transform of the autocorrelation function
                                of a time series.''',
                                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('-f', '--filename', default='data.dat',
                        help='''Either a time series or an acf file.''')

    parser.add_argument('--filetype', default='series',
                        choices=['series', 'acf'],
                        help='''Type of the input file.''')

    parser.add_argument('--c1', default=1, type=int,
                        help='''Time Column.''')

    parser.add_argument('--c2', default=['2'], nargs='+',
                        help='''Data Columns.''')

    parser.add_argument('-tu', '--timeunit',
                        choices=["s", "ms", "mus", "ns", "ps", "fs"],
                        type=str, default="ps",
                        help='''Unit of the time series''')

    parser.add_argument('--show',
                        help='''Show the plot in an external window.''',
                        default=False, action='store_true')

    parser.add_argument('-o', '--output', type=str, default=None,
                        help='''Output File.''')

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

    mean = np.mean(series)
    var = np.var(series)
    ny = series - mean
    acf = np.correlate(ny, ny, "full")[len(ny) - 1 :]
    acf_t = acf / var / len(ny)

    return acf_t


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
    # bl = args.baseline

    data = np.loadtxt(f)

    if args.filetype == 'series':
        x = data[:,0]
        y = data[:,1]
        acf_y = acf(y)

    elif args.filetype == 'acf':
        x = data[:,c1]
        y = None
        acf_y = data[:,c2].flatten()

    #
    # Save to a file
    #
    if args.filetype == 'series':

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

    # This is not needed anymore thanks to the correction of the DCT
    # #
    # # Subtract a baseline from the Cosine Transform before multiplication
    # # with the prefactor (see JPCB, 2013, 117, 7157)
    # #
    # if bl:
    #     idxs = np.where(freqs > bl)
    #     avg = np.mean(specden_ft_part[idxs])
    #     specden_ft_part -= avg
    #     specden_ft_part[specden_ft_part < 0] = 0

    #
    # Calculate the total Spectral Density and convert to wavenumbers
    #
    prefac = freqs / (k_B * T * np.pi)
    specden = prefac * specden_ft_part * eV2wn

    #
    # Save to a file
    #
    data = np.c_[freqs, specden_ft_part,specden]
    if not args.output:
        basename = '.'.join(f.split('.')[:-1])
        outfile = basename + '.specden.out'

    else:
        outfile = args.output + '.specden.out'

    header = "\n Frequency (cm^-1)  FT (eV^2 / cm^-1) Spectral Density (cm^-1)\n"
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
        ax.plot(freqs, specden, label="SpecDen")
        ax.set_xlabel(r'$\omega$ / cm$^{-1}$', size=22)
        ax.set_ylabel(r'$S(\omega)$ / cm$^{-1}$', size=22)

        xtickmaj = ticker.MaxNLocator(5)
        xtickmin = ticker.AutoMinorLocator(5)
        ytickmaj = ticker.MaxNLocator(5)
        ytickmin = ticker.AutoMinorLocator(5)
        ax.xaxis.set_major_locator(xtickmaj)
        ax.xaxis.set_minor_locator(xtickmin)
        ax.yaxis.set_major_locator(ytickmaj)
        ax.yaxis.set_minor_locator(ytickmin)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(axis='both', which='major', direction='in', labelsize=22, pad=10, length=5)
        ax.tick_params(axis='both', which='minor', direction='in', labelsize=22, pad=10, length=2)

        # ax0 = ax.twinx()
        # ax0.plot(freqs, specden_ft_part, color="g", label="FT")
        # ax0.set_ylabel(r'eV$^2$ / cm$^{-1}$', size=22)

        # xtickmaj = ticker.MaxNLocator(5)
        # xtickmin = ticker.AutoMinorLocator(5)
        # ytickmaj = ticker.MaxNLocator(5)
        # ytickmin = ticker.AutoMinorLocator(5)
        # ax0.xaxis.set_major_locator(xtickmaj)
        # ax0.xaxis.set_minor_locator(xtickmin)
        # ax0.yaxis.set_major_locator(ytickmaj)
        # ax0.yaxis.set_minor_locator(ytickmin)
        # ax0.xaxis.set_ticks_position('both')
        # ax0.yaxis.set_ticks_position('both')
        # ax0.tick_params(axis='both', which='major', direction='in', labelsize=22, pad=10, length=5)
        # ax0.tick_params(axis='both', which='minor', direction='in', labelsize=22, pad=10, length=2)

        #
        # Plot ACF
        #
        ax1 = plt.subplot(gs[1])
        ax1.plot(x, acf_y, label="ACF")
        ax1.set_xlabel('$t$ / %s' % args.timeunit, size=22)
        ax1.set_ylabel('$R(t)$ / eV$^2$', size=22)

        xtickmaj = ticker.MaxNLocator(5)
        xtickmin = ticker.AutoMinorLocator(5)
        ytickmaj = ticker.MultipleLocator(0.2)
        ytickmin = ticker.AutoMinorLocator(2)
        ax1.xaxis.set_major_locator(xtickmaj)
        ax1.xaxis.set_minor_locator(xtickmin)
        ax1.yaxis.set_major_locator(ytickmaj)
        ax1.yaxis.set_minor_locator(ytickmin)
        ax1.xaxis.set_ticks_position('both')
        ax1.yaxis.set_ticks_position('both')
        ax1.tick_params(axis='both', which='major', direction='in', labelsize=22, pad=10, length=5)
        ax1.tick_params(axis='both', which='minor', direction='in', labelsize=22, pad=10, length=2)

        if y is not None:

            #
            # Plot time series
            #
            ax2 = plt.subplot(gs[2])
            ax2.plot(x, y * 1e3)
            ax2.set_xlabel('$t$ / %s' % args.timeunit, size=22)
            ax2.set_ylabel('$J$ / meV', size=22)

            xtickmaj = ticker.MaxNLocator(5)
            xtickmin = ticker.AutoMinorLocator(5)
            ytickmaj = ticker.MaxNLocator(5)
            ytickmin = ticker.AutoMinorLocator(5)
            ax2.xaxis.set_major_locator(xtickmaj)
            ax2.xaxis.set_minor_locator(xtickmin)
            ax2.yaxis.set_major_locator(ytickmaj)
            ax2.yaxis.set_minor_locator(ytickmin)
            ax2.xaxis.set_ticks_position('both')
            ax2.yaxis.set_ticks_position('both')
            ax2.tick_params(axis='both', which='major', direction='in', labelsize=22, pad=10, length=5)
            ax2.tick_params(axis='both', which='minor', direction='in', labelsize=22, pad=10, length=2)

        plt.tight_layout()
        plt.show()
