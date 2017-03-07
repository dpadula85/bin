#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd
import argparse as arg
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import norm


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Plots a nice graph for MD properties.', formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('-f', '--filename', default='data.dat', help='''File data.dat from mdanalyzer.''')

    parser.add_argument('--c1', default=1, type=int, help='''Column for the progression of the property.''')

    parser.add_argument('--c2', default=['2'], nargs='+', help='''Columns containing the data.''')

    parser.add_argument('--merge', default=False, action='store_true', help='''Merge the data from all
    specified columns in --c2 option in a cumulative data set.''')

    parser.add_argument('--compare', default=False, action='store_true', help='''Plot all data from all
    specified columns in --c2 option together.''')

    parser.add_argument('-ux', '--unitx', type=str, default=None, help='''Unit of the property to be plotted on X axis.''')

    parser.add_argument('-ty', '--titley', type=str, default=None, help='''Name of the property to be plotted on Y axis.''')

    parser.add_argument('-uy', '--unity', type=str, default=None, help='''Unit of the property to be plotted on Y axis.''')

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


def plot_data(data, yhdrs, xhdr=None, ux=None, ty=None, uy=None):

    # Two subplots, unpack the axes array immediately
    # fig = plt.figure(figsize=(16, 12)) 
    fig = plt.figure(figsize=(11.69, 8.27)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])
    stats = np.array([]).reshape(0,4)

    if xhdr:
        x = data[xhdr]
    else:
        x = np.arange(1, len(data) + 1)

    for yhdr in yhdrs:

        y = data[yhdr]

        avg = np.average(y)
        sigma = np.std(y)
        ymin  = y.min() 
        ymax  = y.max()

        stat = np.array([avg, sigma, ymin, ymax])
        stats = np.vstack((stats, stat))

        # Sturges' formula for number of bins
        nbins = int(np.log2(len(y)) + 1)

        #
        # Trajectory subplot
        #
        ax0 = plt.subplot(gs[0])
        ax0.set_xlabel(xhdr, size=26)
        ax0.set_xlim(x.min(), x.max())
        ax0.tick_params(axis='both', which='major', labelsize=24, pad=10)
        line = ax0.plot(x, y, label=yhdr)
        # line = ax0.plot(x, y, marker='.', label=yhdr)
        clr = line[0].get_color()
        ax0.minorticks_on()

        if xhdr:
            title = xhdr

        else:
            title = "Data Set"
            
        if ux:
            ax0.set_xlabel('%s (%s)' % (title, ux), size=26)
        else:
            ax0.set_xlabel('%s' % title, size=26)

        if ty:
            ty = ty.title()
            if uy:
                ax0.set_ylabel('%s (%s)' % (ty, uy), size=26)
            else:
                ax0.set_ylabel('%s' % ty, size=26)

        # Get y scale to set the same for the histogram
        ylim_low, ylim_high = ax0.get_ylim()

        # Average line
        ax0.axhline(avg, linestyle='dashed', linewidth=2, color=clr)

        #
        # Histogram subplot
        #
        ax1 = plt.subplot(gs[1])
        ax1.set_ylim(ylim_low, ylim_high)
        ax1.set_yticklabels([])
        ax1.minorticks_on()
        ax1.set_xlabel('Count', size=26)
        ax1.tick_params(axis='x', which='major', labelsize=24, pad=10)
        ax1.tick_params(axis='x', which='minor', bottom='off')
        ax1.tick_params(axis='x', which='minor', top='off')
        ax1.axhline(avg, linestyle='dashed', linewidth=2, color=clr)

        # Distribution histograms, the graph will be rotated by 90 deg
        n, bins, patches = ax1.hist(y, bins=nbins, orientation='horizontal',
                                    histtype='bar', rwidth=0.75, hatch='//',
                                    fill=False, color=clr, edgecolor=clr)

        # Fit a gaussian, scaled to the real distribution of the data and add it to the legend
        scale_factor = (bins[1] - bins[0]) * len(y)
        fitx = np.linspace(ymin, ymax, 1000)
        gau_fit = norm.pdf(fitx, avg, sigma) * scale_factor
        gau_line = ax1.plot(gau_fit, fitx, '-', linewidth=2, color=clr)


    ax0.legend(bbox_to_anchor=(0.75, 1.06), loc=10, ncol=len(yhdrs),
               borderaxespad=0, fontsize=24).draw_frame(False)

    # plt.tight_layout()

    return fig, stats


if __name__ == '__main__':

    args = options()

    f = args.filename
    ux = args.unitx
    ty = args.titley
    uy = args.unity
    basename = '.'.join(f.split('.')[:-1])
    checkfile(f)

    # Get columns to process and convert it to python numeration
    c1 = args.c1 - 1
    c2 = map(lambda x: x - 1, extend_compact_list(args.c2))

    data = gen_pandas_df(f)

    if len(data.columns) > 1:
        xhdr = data.columns[c1]
        yhdrs = [ data.columns[i] for i in c2 ]

    else:
        xhdr = None
        yhdrs = data.columns

    if args.compare:

        fig, stat = plot_data(data, yhdrs, xhdr, ux, ty, uy)

        # Save plot as vector image
        if args.save:
        
            print(" > Saving plot for COLS %d-%d..." % (min(c2) + 1, max(c2) + 1))
            print
            plt.savefig('%s_cols%d-%d.%s' % (basename, min(c2) + 1, max(c2) + 1, args.save),
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
        
            print(" > Showing plot for COLS %d-%d..." % (min(c2) + 1, max(c2) + 1))
            print
            plt.show()

    
    if args.merge:

        tot = pd.concat([ data[yhdr] for yhdr in yhdrs ], ignore_index=True).to_frame(name="Merged")

        yhdrs = ["Merged"]
        xhdrs = None
        fig, stat = plot_data(tot, yhdrs, xhdrs, ux, ty, uy)

        # Save plot as vector image
        if args.save:
        
            print(" > Saving plot for %s..." % yhdrs[0])
            print
            plt.savefig('%s_%s.%s' % (basename, yhdrs[0], args.save),
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
        
            print(" > Showing plot for %s..." % yhdrs[0])
            print
            plt.show()


    if not args.merge and not args.compare:

        stat = np.array([]).reshape(0,4)
        for yhdr in yhdrs:

            fig, statcol = plot_data(data, [yhdr], xhdr, ux, ty, uy)
            stat = np.vstack((stat, statcol))

            # Save plot as vector image
            if args.save:

                print(" > Saving plot for %s..." % yhdr)
                print
                plt.savefig('%s_%s.%s' % (basename, yhdr, args.save),
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

                print(" > Showing plot for %s..." % yhdr)
                print
                plt.show()

    print
    print(banner(ch="=", length=60))
    print("Statistical Analysis of %s" % args.filename)
    print
    print("                      Avg.  Std. Dev.       Min.       Max.")
    print(banner(ch="-", length=60))

    for i, yhdr in enumerate(yhdrs):

        data = [ yhdr, stat[i,0], stat[i,1], stat[i,2], stat[i,3]]
        print("%-15s %10.4f %10.4f %10.4f %10.4f" % tuple(data))

    print(banner(ch="=", length=60))
    print
