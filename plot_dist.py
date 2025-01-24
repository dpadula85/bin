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

    inp.add_argument('--mm', type=str, dest='MMFile',
                     help='''MM Data.''')

    inp.add_argument('-c', '--color', type=str, default="b", dest='Color',
                     help='''Color of the plot.''')

    inp.add_argument('-l', '--label', type=str, default="phi_1", dest='Label',
                     help='''Label of the plot.''')


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

def myround(x, base=5):
    return base * round(x/base)

def plot(x, **kwargs):

    color = kwargs.pop("Color", "b")
    label = kwargs.pop("Label", "phi_1")
    nbins = int(np.log2(x.shape[0]) + 1)

    fig = plt.figure(figsize=(12, 9))
    gs = gridspec.GridSpec(1, 1)

    # Plot
    ax = plt.subplot(gs[0])
    n, bins, patches = ax.hist(x, bins=nbins, histtype='bar', rwidth=0.75,
                               hatch='//', fill=False, color=color,
                               edgecolor=color)

    # Fancy stuff
    ax.set_xlabel(r"$\%s_{%s}$ / deg" % tuple(label.split("_")), size=42, labelpad=5)
    ax.set_ylabel(u"Count", size=42)

    xtickmaj = ticker.MultipleLocator(60)
    xtickmin = ticker.AutoMinorLocator(4)
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
    # xmin = myround(x.min(), base=30)
    # xmax = myround(x.max(), base=30)
    ax.set_xlim(-180, 180)

    return

if __name__ == '__main__':

    Opts = options()
    mmdata = np.loadtxt(Opts["MMFile"])[:,:2]
    plot(mmdata[:,1], color=Opts["Color"], label=Opts["Label"])
    plot(mmdata[:,1], **Opts)

    if Opts["Out"]:
        if Opts["OutPre"]:
            name = "%s_dist_%s_%s.%s" % (Opts["OutPre"],
                                         Opts["Label"].split("_")[0],
                                         Opts["Label"].split("_")[1],
                                         Opts["Out"])
        else:
            name = "dist_%s_%s.%s" % (Opts["Label"].split("_")[0],
                                      Opts["Label"].split("_")[1],
                                      Opts["Out"])

        plt.subplots_adjust(bottom=0.15, top=0.9, left=0.175)
        plt.savefig(name, dpi=600)
    else:
        plt.show()
