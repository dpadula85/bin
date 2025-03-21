#!/usr/bin/env python

import sys
import numpy as np
import argparse as arg
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
    # "pgf.preamble": [
    #         r"\usepackage[utf8x]{inputenc}",
    #         r"\usepackage[T1]{fontenc}",
    #         r"\usepackage{siunitx}",
    #     ]
    }

mpl.rcParams.update(pgf_with_latex)



def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument('-f', '--file', type=str, dest='DuschFile',
                     help='''Duschinsky matrix file''')

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


def plot_dusch(J):

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[ 30, 1 ])

    colormap = plt.cm.bwr_r
    normalize = mpl.colors.Normalize(vmin=-1, vmax=1)

    # Plot
    ax = plt.subplot(gs[0])

    mat = ax.matshow(J, cmap=colormap, norm=normalize)

    ax.set_xlabel(r"Normal Mode", size=32, labelpad=5)
    ax.set_ylabel(r"Normal Mode", size=32, labelpad=5)

    xtickmaj = ticker.MultipleLocator(20)
    xtickmin = ticker.AutoMinorLocator(5)
    ax.xaxis.set_major_locator(xtickmaj)
    ax.xaxis.set_minor_locator(xtickmin)
    ax.yaxis.set_major_locator(xtickmaj)
    ax.yaxis.set_minor_locator(xtickmin)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='major', direction='in', labelsize=28, pad=10, length=5, labelbottom=True, labeltop=False)
    ax.tick_params(axis='both', which='minor', direction='in', labelsize=28, pad=10, length=2)

    ax1 = plt.subplot(gs[1])
    cb = plt.colorbar(mat, cax=ax1)
    xtickmaj = ticker.MultipleLocator(0.25)
    ax1.yaxis.set_major_locator(xtickmaj)
    ax1.yaxis.set_ticks_position('both')
    ax1.tick_params(axis='y', which='major', direction='in', labelsize=28, pad=10, length=5, labelright=True, labelleft=False)
    ax1.set_ylim(-1, 1)

    return


if __name__ == '__main__':

    Opts = options()
    J = np.loadtxt(Opts["DuschFile"]).reshape(156,156)
    plot_dusch(J)

    if Opts["Out"]:
        if Opts["OutPre"]:
            name = "%s_dusch.%s" % (Opts["OutPre"], Opts["Out"])
        else:
            name = "dusch.%s" % Opts["Out"]
        plt.subplots_adjust(bottom=0.2, top=0.9)
        plt.savefig(name, dpi=600)
    else:
        plt.show()

    # plot_overlap(ovlp, color=Opts["Color"])

    # if Opts["Out"]:
    #     if Opts["OutPre"]:
    #         name = "%s_ovlp.%s" % (Opts["OutPre"], Opts["Out"])
    #     else:
    #         name = "ovlp.%s" % Opts["Out"]
    #     plt.subplots_adjust(bottom=0.1, top=0.8)
    #     plt.savefig(name, dpi=600)
    # else:
    #     plt.show()
