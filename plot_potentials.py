#!/usr/bin/env python

import sys
import numpy as np
import argparse as arg
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker, gridspec
from matplotlib import rc

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

    inp.add_argument('--qm', type=str, dest='QMFile',
                     help='''QM Data.''')

    inp.add_argument('--mm', type=str, dest='MMFile',
                     help='''MM Data.''')

    inp.add_argument('-m', '--min', type=int, default=180, dest='Min',
                     help='''Centre of the plot.''')

    inp.add_argument('-l', '--label', type=str, default="phi_1", dest='Label',
                     help='''Label of the plot.''')

    inp.add_argument('-c', '--color', type=str, default="b", dest='Color',
                     help='''Color of the plot.''')

    inp.add_argument('-r', '--rescale', default=False, action="store_true", dest='Rescale',
                     help='''Vertical shifting activated.''')

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


def myround(x, base=30):
    return base * round(x/base)


def plot(qmdata, mmdata, **kwargs):

    color = kwargs.pop("Color", "b")
    label = kwargs.pop("Label", "phi_1")
    rescale = kwargs.pop("Rescale", False)
    x_min = kwargs.pop("Min", 180)
    x_min = myround(x_min)

    fig = plt.figure(figsize=(12, 9))
    gs = gridspec.GridSpec(1, 1)

    # Sort and transform
    qmdata = qmdata[qmdata[:,0].argsort()]
    mmdata = mmdata[mmdata[:,0].argsort()]
    mask = qmdata[:,0] < 0
    qmdata[:,0][mask] = 360 + qmdata[:,0][mask]
    mask = mmdata[:,0] < 0
    mmdata[:,0][mask] = 360 + mmdata[:,0][mask]
    qmdata = qmdata[qmdata[:,0].argsort()]
    mmdata = mmdata[mmdata[:,0].argsort()]
    # qmdata[:,0] = np.unwrap(qmdata[:,0] * np.pi / 180) * 180 / np.pi - 360
    # mmdata[:,0] = np.unwrap(mmdata[:,0] * np.pi / 180) * 180 / np.pi - 360

    # Vertical shift
    if rescale:
        qmdata[:,1] = qmdata[:,1] - qmdata[:,1].min()
        mmdata[:,1] = mmdata[:,1] - mmdata[:,1].min()

    # Plot
    ax = plt.subplot(gs[0])

    ax.scatter(qmdata[:,0], qmdata[:,1], color=color, s=80)
    ax.plot(mmdata[:,0], mmdata[:,1], color=color, lw=2, ls="--")

    ax.set_xlabel(r"$\%s_{%s}$ / deg" % tuple(label.split("_")), size=42, labelpad=5)
    ax.set_ylabel(r"$\Delta E$ / KJ mol$^{-1}$", size=42, labelpad=5)

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
    ax.tick_params(axis='both', which='major', direction='in', labelsize=38, pad=10, length=5)
    ax.tick_params(axis='both', which='minor', direction='in', labelsize=38, pad=10, length=2)
    ax.set_xlim(qmdata[:,0].min(), qmdata[:,0].max())
    ax.set_xlim(x_min - 90, x_min + 90)

    return


if __name__ == '__main__':

    Opts = options()
    qmdata = np.loadtxt(Opts["QMFile"])[:,:2]
    mmdata = np.loadtxt(Opts["MMFile"])[:,:2]
    plot(qmdata, mmdata, **Opts)

    if Opts["Out"]:

        if Opts["OutPre"]:
            name = "%s_%s_%s.%s" % (Opts["OutPre"],
                                    Opts["Label"].split("_")[0],
                                    Opts["Label"].split("_")[1],
                                    Opts["Out"])
        else:
            name = "%s_%s.%s" % (Opts["Label"].split("_")[0],
                                 Opts["Label"].split("_")[1],
                                 Opts["Out"])

        plt.subplots_adjust(bottom=0.15, top=0.9, left=0.15)
        plt.savefig(name, dpi=600)
    else:
        plt.show()
