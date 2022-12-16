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

    inp.add_argument('-f', '--file', type=str, dest='ModesFile',
                     help='''Normal Modes comparison file.''')

    inp.add_argument('-c', '--color', type=str, default="b", dest='Color',
                     help='''Color of the plot.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-p', '--pre', default=None, type=str, dest='OutPre',
                     help='''Output File Prefix.''')

    out.add_argument('-o', '--output', default=None, type=str, dest='Out',
                     help='''Output File Format.''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


def plot_modes(modes, color="b"):

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 1)

    # Plot
    ax = plt.subplot(gs[0])
    x = np.arange(0, modes.max() + modes.max() / 10)
    y = x

    ax.plot(x, y, color="k", ls="--")
    ax.scatter(np.sort(modes[:,0]), np.sort(modes[:,1]), color=color)

    ax.set_xlabel(r"$\nu^{QM}$ / cm$^{-1}$", size=18, labelpad=5)
    ax.set_ylabel(r"$\nu^{MM}$ / cm$^{-1}$", size=18, labelpad=5)

    xtickmaj = ticker.MultipleLocator(500)
    xtickmin = ticker.AutoMinorLocator(5)
    ytickmaj = ticker.MultipleLocator(500)
    ytickmin = ticker.AutoMinorLocator(5)
    ax.xaxis.set_major_locator(xtickmaj)
    ax.xaxis.set_minor_locator(xtickmin)
    ax.yaxis.set_major_locator(ytickmaj)
    ax.yaxis.set_minor_locator(ytickmin)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='major', direction='in', labelsize=16, pad=10, length=5)
    ax.tick_params(axis='both', which='minor', direction='in', labelsize=16, pad=10, length=2)
    ax.set_xlim(0, modes.max() + 50)
    ax.set_ylim(0, modes.max() + 50)

    return


def plot_overlap(overlap, color="b"):

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 1)

    # Sort
    overlap = overlap[overlap[:,0].argsort()]

    # Plot
    ax = plt.subplot(gs[0])

    ax.bar(overlap[:,0], overlap[:,1], color=color, width=1.0)

    ax.set_xlabel(r"Normal Mode", size=18, labelpad=5)
    # ax.set_ylabel(r"overlap", size=18, labelpad=5)
    ax.set_ylabel(r"$\langle Q_i^{MM} \vert Q_i^{QM} \rangle$", size=18, labelpad=5)
    ax.xaxis.set_label_position('top')

    xtickmaj = ticker.MultipleLocator(20)
    xtickmin = ticker.AutoMinorLocator(5)
    ytickmaj = ticker.MultipleLocator(0.2)
    ytickmin = ticker.AutoMinorLocator(3)
    ax.xaxis.set_major_locator(xtickmaj)
    ax.xaxis.set_minor_locator(xtickmin)
    ax.yaxis.set_major_locator(ytickmaj)
    ax.yaxis.set_minor_locator(ytickmin)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='major', direction='in', labelsize=16, pad=10, length=5, labeltop=True, labelbottom=False)
    ax.tick_params(axis='both', which='minor', direction='in', labelsize=16, pad=10, length=2)
    ax.set_xlim(0, overlap[:,0].max())
    ax.set_ylim(0, 1)

    return


if __name__ == '__main__':

    Opts = options()
    modes_file = np.loadtxt(Opts["ModesFile"])
    modes = np.c_ [ modes_file[:,3], modes_file[:,1] ]
    ovlp = np.c_ [ modes_file[:,0], modes_file[:,-2] ]
    plot_modes(modes, color=Opts["Color"])

    if Opts["Out"]:
        if Opts["OutPre"]:
            name = "%s_modes.%s" % (Opts["OutPre"], Opts["Out"])
        else:
            name = "modes.%s" % Opts["Out"]
        plt.savefig(name, dpi=600)
    else:
        plt.show()

    plot_overlap(ovlp, color=Opts["Color"])

    if Opts["Out"]:
        if Opts["OutPre"]:
            name = "%s_ovlp.%s" % (Opts["OutPre"], Opts["Out"])
        else:
            name = "ovlp.%s" % Opts["Out"]
        plt.savefig(name, dpi=600)
    else:
        plt.show()
