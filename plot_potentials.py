#!/usr/bin/env python

import sys
import numpy as np
import argparse as arg
import matplotlib.pyplot as plt
from matplotlib import ticker, gridspec
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument('--qm', type=str, dest='QMFile',
                     help='''QM Data.''')

    inp.add_argument('--mm', type=str, dest='MMFile',
                     help='''MM Data.''')

    inp.add_argument('-l', '--label', type=int, default=1, dest='Label',
                     help='''Label of the plot.''')

    inp.add_argument('-c', '--color', type=str, default="b", dest='Color',
                     help='''Color of the plot.''')

    inp.add_argument('-r', '--rescale', default=False, action="store_true", dest='Rescale',
                     help='''Vertical shifting activated.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--output', default=None, type=str, dest='Out',
                     help='''Output File Format.''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


def plot(qmdata, mmdata, color="b", label=1, rescale=False):

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(1, 1)

    # Sort and transform
    qmdata = qmdata[qmdata[:,0].argsort()]
    mmdata = mmdata[mmdata[:,0].argsort()]
    mask = qmdata[:,0] < 0
    qmdata[:,0][mask] = 360 + qmdata[:,0][mask]

    # Vertical shift
    if rescale:
        qmdata[:,1] = qmdata[:,1] - qmdata[:,1].min()
        mmdata[:,1] = mmdata[:,1] - mmdata[:,1].min()

    # Plot
    ax = plt.subplot(gs[0])

    ax.scatter(qmdata[:,0], qmdata[:,1], color=color)
    ax.plot(mmdata[:,0], mmdata[:,1], color=color, lw=1, ls="--")

    ax.set_xlabel(r"$\delta_{%d}$ / deg" % label, size=18, labelpad=5)
    ax.set_ylabel(r"$\Delta E$ / KJ mol$^{-1}$", size=18, labelpad=5)

    xtickmaj = ticker.MultipleLocator(30)
    xtickmin = ticker.AutoMinorLocator(5)
    ytickmaj = ticker.MaxNLocator(5)
    ytickmin = ticker.AutoMinorLocator(5)
    ax.xaxis.set_major_locator(xtickmaj)
    ax.xaxis.set_minor_locator(xtickmin)
    ax.yaxis.set_major_locator(ytickmaj)
    ax.yaxis.set_minor_locator(ytickmin)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='major', direction='in', labelsize=16, pad=10, length=5)
    ax.tick_params(axis='both', which='minor', direction='in', labelsize=16, pad=10, length=2)
    ax.set_xlim(qmdata[:,0].min(), qmdata[:,0].max())

    return


if __name__ == '__main__':

    Opts = options()
    qmdata = np.loadtxt(Opts["QMFile"])[:,:2]
    mmdata = np.loadtxt(Opts["MMFile"])[:,:2]
    plot(qmdata, mmdata,
         color=Opts["Color"],
         label=Opts["Label"],
         rescale=Opts["Rescale"])

    if Opts["Out"]:
        plt.savefig("delta_%d.%s" % (Opts["Label"], Opts["Out"]), dpi=600)
    else:
        plt.show()
