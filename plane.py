#!/usr/bin/env python

import sys
import textwrap
import numpy as np
import argparse as arg
from os.path import splitext
from scipy.optimize import minimize

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Fits a plane to a molecule',
            formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('-f', '--filename', type=str, required=True,
            help='''File containing the molecular geometry in G09 format.''')

    parser.add_argument('-s', '--sel', nargs='+', type=str, default=None,
            help='''Atom selection for the plane calculation.''')

    args = parser.parse_args()

    return args


def process_selection(string):

    string =  ','.join(string).replace(',,',',')

    try:
        f = open(string, 'r')
        string = f.readlines()
        f.close()
        string =  ','.join(string).replace(',,',',')
        string = string.replace(',', ' ')
        string = map(lambda x: x - 1, extend_compact_list(string))

    except IOError:
        string = string.replace(',', ' ')
        string = map(lambda x: x - 1, extend_compact_list(string))

    return string


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    idxs = idxs.split()

    for idx in idxs:

        to_extend = idx.split('-')

        if len(to_extend) > 1:

            sel =  map(int, to_extend)
            extended += range(sel[0],sel[1]+1,1)

        else:
        
            extended.append(int(idx))
    
    return extended


def format_selection(intlist):

    s = ''
    for i in intlist:
        s += '%2d ' % (i + 1) 

    return s


def parse_log(logfile):

    structure = []
    opt = 0
    with open(logfile, 'r') as f:
        for line in f:
            if "Optimization completed" in line:
                opt = 1
    
            # It could be Standard or Input orientation
            # depending on nosymm keyword
            if "orientation:" in line and opt == 1:
                # Skip the head-of-table lines (4 lines)
                next(f)
                next(f)
                next(f)
                next(f)
                curr_line = next(f).split()
                while len(curr_line) == 6:
                    atom_n = int(curr_line[1])
                    atom_x = float(curr_line[3])
                    atom_y = float(curr_line[4])
                    atom_z = float(curr_line[5])
                    structure.append([atom_n, atom_x, atom_y, atom_z])
                    curr_line = next(f).split()
                opt = 0
    
    return np.array(structure)


def plane(x, y, parms):

    a = parms[0]
    b = parms[1]
    c = parms[2]
    z = a*x + b*y + c

    return z


def error(parms, x, y, z):

    target = plane(x, y, parms)
    err = np.sqrt(np.mean((target - z)**2))

    return err


if __name__ == "__main__":

    args = options()


    # Read Data
    outpref = splitext(args.filename)[0]

    struct = parse_log(args.filename)
    with open("%s.xyz" % outpref, 'w') as f:
        f.write("%d\n" % len(struct))
        f.write("\n")
        np.savetxt(f, struct, fmt="%3d %14.8f %14.8f %14.8f")
    
    if args.sel:
        idxs = process_selection(args.sel)
        points = struct[idxs,1:]
    else:
        points = struct[:,1:]

    x = points[:,0]
    y = points[:,1]
    z = points[:,2]

    # Calculate the Plane
    guess = np.zeros(3)
    res = minimize(error, guess, args=(x, y, z))
    
    a = res.x[0]
    b = res.x[1]
    c = res.x[2]
    err = res.fun

    print
    print("--------- Results of the Fitting Procedure ---------") 
    print
    print("Input: %s" % args.filename)

    if args.sel:
        print("Selection:\n%s" % textwrap.fill(format_selection(idxs), width=55))

    print
    print("RMSD: %13.7f" % err)
    print("Equation: z = %13.7f*x + %13.7f*y + %13.7f" % (a, b, c)) 
    print
    print("Saving the vmd script %s.vmd" % outpref)
    print
    
    # Generate a plottable surface
    point  = np.array([0.0, 0.0, c])
    normal = np.array(np.cross([1,0,a], [0,1,b]))
    d = -point.dot(normal)
    xx, yy = np.meshgrid([x.min(), x.max()], [y.min(), y.max()])
    zz = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

    # Save vmd script plotting the plane and the molecule
    with open("%s.vmd" % outpref, 'w') as f:

        f.write("mol new %s.xyz type xyz\n" % outpref)
        f.write("draw color yellow\n")
        f.write("draw triangle {%7.4f %7.4f %7.4f} {%7.4f %7.4f %7.4f} {%7.4f %7.4f %7.4f}\n" %
                (xx[0][0], yy[0][0], zz[0][0], xx[0][1], yy[0][1], zz[0][1], xx[1][0], yy[1][0], zz[1][0]))
        f.write("draw triangle {%7.4f %7.4f %7.4f} {%7.4f %7.4f %7.4f} {%7.4f %7.4f %7.4f}\n" %
                (xx[1][1], yy[1][1], zz[1][1], xx[1][0], yy[1][0], zz[1][0], xx[0][1], yy[0][1], zz[0][1]))
