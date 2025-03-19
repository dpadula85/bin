#!/usr/bin/env python

import sys
import numpy as np
import argparse as arg
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker, gridspec
from matplotlib import rc

from scipy.constants import k as k_b
from scipy.constants import N_A

k_b = k_b * 1e-3 * N_A # in kj / mol

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

    inp.add_argument('-n', '--nats', type=int, dest='NAts',
                     help='''Number of atoms for virial computation.''')

    inp.add_argument('-t', '--temp', type=int, dest='T',
                     help='''Temperature for virial computation.''')

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
    T = kwargs.pop("T", 298.0)
    NAts = kwargs.pop("NAts", 0)
    x_min = myround(x_min)

    virial = (3 * NAts - 6) * k_b * T * 0.5

    fig = plt.figure(figsize=(12, 9))
    gs = gridspec.GridSpec(1, 1)

    # Sort and transform
    qmdata = qmdata[qmdata[:,0].argsort()][1:]
    mmdata = mmdata[mmdata[:,0].argsort()][1:]

    # Plot
    ax = plt.subplot(gs[0])

    ax.axhline(virial, color="r", lw=2, ls="--", label=r"$\frac{(3N - 6)}{2} \cdot k_B T$")
    ax.plot(qmdata[:,0], qmdata[:,1], color="g", marker="o", markersize=8, label=r"QM @ QMD--FF")
    ax.plot(mmdata[:,0], mmdata[:,1], color="orange", marker="o", markersize=8, label=r"QM @ GAFF")

    ax.set_xlabel(r"$t$ / ps", size=42, labelpad=5)
    ax.set_ylabel(r"$\Delta E$ / kJ mol$^{-1}$", size=42, labelpad=5)

    xtickmaj = ticker.MultipleLocator(100)
    xtickmin = ticker.AutoMinorLocator(2)
    ytickmaj = ticker.MultipleLocator(100)
    ytickmin = ticker.AutoMinorLocator(5)
    ax.xaxis.set_major_locator(xtickmaj)
    ax.xaxis.set_minor_locator(xtickmin)
    ax.yaxis.set_major_locator(ytickmaj)
    ax.yaxis.set_minor_locator(ytickmin)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='major', direction='in', labelsize=38, pad=10, length=5)
    ax.tick_params(axis='both', which='minor', direction='in', labelsize=38, pad=10, length=2)
    ax.set_ylim(50, 550)
    ax.set_xlim(100, 500)
    ax.legend(fontsize=32, frameon=False, loc="upper left")

    return


if __name__ == '__main__':

    Opts = options()
    qmdata = np.loadtxt(Opts["QMFile"])[:,:2]
    mmdata = np.loadtxt(Opts["MMFile"])[:,:2]
    plot(qmdata, mmdata, **Opts)

    if Opts["Out"]:

        if Opts["OutPre"]:
            name = "%s_virial.%s" % (Opts["OutPre"],
                                     Opts["Out"])
        else:
            name = "virial.%s" % (Opts["Out"])

        plt.subplots_adjust(bottom=0.15, top=0.9, left=0.15)
        plt.savefig(name, dpi=600)
    else:
        plt.show()
