#!/usr/bin/env python

import os
import sys
import numpy as np
import argparse as arg
import matplotlib.pyplot as plt

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Plots a nice graph with excitonic coefficients.')

    # Optional arguments
    parser.add_argument('-f', '--filename', default='diag.dat', help='''File diag.dat from EXAT.''')

    parser.add_argument('-s', '--save', help='''Save the plot as an image.
    Specify the extension.''')

    parser.add_argument('--show', help='''Show the plot in an external window.''',
    default=False, action='store_true')

    parser.add_argument('--square', help='''Use squares of coefficients for the plot.''',
    default=False, action='store_true')

    args = parser.parse_args()

    if not args.show and not args.save:
        parser.print_help()
        sys.exit()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return args


def doplot(coeff):

    fig, ax = plt.subplots()
    im = ax.pcolor(coeff)
    fig.colorbar(im)
    
    # Generate labels
    xlabels = np.arange(1, len(coeff[0])+1)
    
    # X axis
    # Hide main labels and assign to minor labels their value
    plt.xlabel("Local States")
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xlabels, visible=False)
    ax.set_xticks(xlabels - 0.5, minor=True)
    ax.set_xticklabels(xlabels, minor=True)
    ax.set_xlim(right=xlabels[-1])
    
    # Y axis
    # Hide main labels and assign to minor labels their value
    plt.ylabel("Excitonic States")
    ax.set_yticks(xlabels)
    ax.set_yticklabels(xlabels, visible=False)
    ax.set_yticks(xlabels - 0.5, minor=True)
    ax.set_yticklabels(xlabels, minor=True)
    ax.set_ylim(top=xlabels[-1])
    
    # Rotate labels for a fancier plot
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=-45)
    plt.setp(ax.yaxis.get_minorticklabels(), rotation=-45)
    
    # Show grid for separation of matrix elements
    plt.grid()


def checkfile(filename):

    if not os.path.isfile(filename):
        print("\nFile %s not found!\n" % filename)
        sys.exit()


if __name__ == '__main__':

    args = options()

    if args.save:
        imgfmt = args.save
    else:
        imgfmt = False

    f = args.filename
    checkfile(f)
    diag = np.loadtxt(f)
    
    # Get coefficients from the diagonalized matrix
    coeff = diag[:diag.shape[0]/2,2:]
    
    if args.square:
        coeff = np.square(coeff)

    doplot(coeff)

    # Save plot as vector image
    if imgfmt:
        plt.savefig('coeffs.%s' % imgfmt)
    
    # Show the plot
    if args.show:
        plt.show()
