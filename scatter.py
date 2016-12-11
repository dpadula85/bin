#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd
import argparse as arg
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib import gridspec


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Plots a series of nice scatter plots.', formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('-f', '--filename', default='data.dat', help='''File data.dat from mdanalyzer.''')

    parser.add_argument('-c', '--cols', default=['1', '2'], nargs='+', help='''Columns containing the data.''')

    parser.add_argument('-t', '--title', type=str, default=None, help='''Name of the property to be plotted on Y axis.''')

    parser.add_argument('-u', '--unit', type=str, default=None, help='''Unit of the property to be plotted on Y axis.''')

    parser.add_argument('-s', '--save', help='''Save the plot as an image. Specify the extension.''')

    parser.add_argument('--show', help='''Show the plot in an external window.''',
    default=False, action='store_true')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

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


def checkfile(filename):

    if not os.path.isfile(filename):
        print(banner(text='ERROR', ch='=', length=80))
        print(" File %s not found!" % filename)
        sys.exit()


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


def gen_pandas_df(filename):

    df = pd.read_csv(filename, delim_whitespace=True, comment='#', header=None)

    with open(filename) as f:
        for line in f:
            if line.startswith('#') and len(line.strip()) > 1:
                header = line.split()[1:]

    df.columns = header

    return df


def plot_data(data, hdrs, t=None, u=None):

    # Two subplots, unpack the axes array immediately
    # fig = plt.figure(figsize=(16, 12)) 
    nplots = (len(hdrs) * len(hdrs) - 1) / 2
    fig = plt.figure(figsize=(11.69, 8.27))

    if nplots > 1:
        gs = gridspec.GridSpec(nplots / 3, 3)

    else:
        gs = gridspec.GridSpec(1, 1)

    k = 0
    for i in range(len(hdrs)):
        for j in range(i, len(hdrs)):

            if i != j:

                x = data[hdrs[i]]
                y = data[hdrs[j]]

                ax = plt.subplot(gs[k])
                ax.scatter(x, y, marker='o', color='b')

                lims = [np.min([ax.get_xlim(), ax.get_ylim()]),
                        np.max([ax.get_xlim(), ax.get_ylim()])]

                # Best Fitting Line
                p, V = np.polyfit(x, y, 1, cov=True)
                fit_x = np.linspace(np.min(ax.get_xlim()), np.max(ax.get_xlim()), 100)
                fit_y = p[0] * fit_x + p[1]
                xerr = np.sqrt(V[0,0])
                yerr = np.sqrt(V[1,1])
                r_squared = 1 - (sum((y - (p[0] * x + p[1]))**2) / ((len(y) - 1) * np.var(y, ddof=1)))
                ax.plot(fit_x, fit_y, color='b', lw=1.5, ls='--')

                print
                print(banner(ch="=", length=60))
                print("Line Fitting of %s vs %s" % (hdrs[i], hdrs[j]))
                print
                print(banner(ch="-", length=60))
                print("             Value            Error")
                print(banner(ch="-", length=60))
                print("slope %16.6e %16.6e" % (p[0], xerr))
                print("intcp %16.6e %16.6e" % (p[1], yerr))
                print("R^2   %12.6f" % r_squared)
                print(banner(ch="=", length=60))

                # x = y line
                ax.plot(lims, lims, color='k')

                ax.set_aspect('equal')
                ax.set_xlim(lims)
                ax.set_ylim(lims)
                tickmaj = ticker.MaxNLocator(5)
                tickmin = ticker.AutoMinorLocator(5)
                ax.xaxis.set_major_locator(tickmaj)
                ax.xaxis.set_minor_locator(tickmin)
                ax.yaxis.set_major_locator(tickmaj)
                ax.yaxis.set_minor_locator(tickmin)

                ax.tick_params(axis='both', which='major', labelsize=16, pad=10)
                ax.minorticks_on()

                if t:
                    if u:
                        ax.set_xlabel('%s\n%s(%s)' % (hdrs[i], t, u), size=22)
                        ax.set_ylabel('%s\n%s(%s)' % (hdrs[j], t, u), size=22)

                    else:
                        ax.set_xlabel('%s\n%s' % (hdrs[i], t), size=22)
                        ax.set_ylabel('%s\n%s' % (hdrs[j], t), size=22)

                else:
                    ax.set_xlabel('%s' % hdrs[i], size=22)
                    ax.set_ylabel('%s' % hdrs[j], size=22)

                k += 1
                
    # ax0.legend(bbox_to_anchor=(0.75, 1.06), loc=10, ncol=len(yhdrs),
    #            borderaxespad=0, fontsize=24).draw_frame(False)

    plt.tight_layout()

    return fig


if __name__ == '__main__':

    args = options()

    f = args.filename
    t = args.title
    u = args.unit
    basename = '.'.join(f.split('.')[:-1])
    checkfile(f)

    # Get columns to process and convert it to python numeration
    cols = map(lambda x: x - 1, extend_compact_list(args.cols))

    data = gen_pandas_df(f)
    hdrs = [ data.columns[i] for i in cols ]

    fig = plot_data(data, hdrs, t, u)

    if args.save:
    
        print(" > Saving plot ...")
        print
        plt.savefig('%s_scatter.%s' % (basename, args.save),
                    dpi=1200, transparent=True, bbox_inches='tight')
    
    # Show the plot
    if args.show:
    
        # Uncomment the two linex of the backend in use to generate a
        # maximized-window plot
    
        # Option 1
        # QT backend
        # manager = plt.get_current_fig_manager()
        # manager.window.showMaximized()
    
        # Option 2
        # TkAgg backend
        # manager = plt.get_current_fig_manager()
        # manager.resize(*manager.window.maxsize())
    
        # Option 3
        # WX backend
        # manager = plt.get_current_fig_manager()
        # manager.frame.Maximize(True)
    
        print(" > Showing plot ...")
        print
        plt.show()
