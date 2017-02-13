#!/usr/bin/env python

import sys
import time
import numpy as np
import argparse as arg
from elements import ELEMENTS
from functools import partial
import multiprocessing as mp

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Calculates the volume of a molecule with a Monte Carlo approach.', formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('-f', '--filename', required=True, help='''XYZ structure file.''')

    parser.add_argument('-n', '--npts', type=int, default=1e5, help='''Number of points for Monte Carlo sampling.''')

    args = parser.parse_args()

    return args


def make_box(coords, VdWs):

    xmin = 0.0
    xmax = 0.0
    ymin = 0.0
    ymax = 0.0
    zmin = 0.0
    zmax = 0.0

    for i, coord in enumerate(coords):

        xm = coord[0] - VdWs[i]
        xM = coord[0] + VdWs[i]
        ym = coord[1] - VdWs[i]
        yM = coord[1] + VdWs[i]
        zm = coord[2] - VdWs[i]
        zM = coord[2] + VdWs[i]

        if xm < xmin:
            xmin = xm

        if xM > xmax:
            xmax = xM

        if ym < ymin:
            ymin = ym

        if yM > ymax:
            ymax = yM

        if zm < zmin:
            zmin = zm

        if zM > zmax:
            zmax = zM

    box_edges = np.array([[xmin, xmax],
                          [ymin, ymax],
                          [zmin, zmax]])

    return box_edges


def inside_VdW(coords, VdWs, point):

    for i, coord in enumerate(coords):
        r = VdWs[i]
        d = np.linalg.norm(point - coord)
        
        if d <= r:
            return True

    return False


def calc_vol(coords, VdWs, npts=100000):

    box = make_box(coords, VdWs)
    box_x = box[0,:]
    box_y = box[1,:]
    box_z = box[2,:]

    lx = box_x[1] - box_x[0]
    ly = box_y[1] - box_y[0]
    lz = box_z[1] - box_z[0]
    box_vol = lx * ly * lz

    p = mp.Pool(mp.cpu_count())
    iterable = range(int(npts))
    func = partial(monte_carlo_hit, box, coords, VdWs)
    hits = p.map(func, iterable)
    pts_in = sum(hits)

    vol_ratio = float(pts_in) / npts

    return vol_ratio * box_vol


def monte_carlo_hit(box, coords, VdWs, seed=None):

    np.random.seed(seed + int(time.time()))

    box_x = box[0,:]
    box_y = box[1,:]
    box_z = box[2,:]

    x = np.random.uniform(box_x[0], box_x[1])
    y = np.random.uniform(box_y[0], box_y[1])
    z = np.random.uniform(box_z[0], box_z[1])
    pt = np.array([x, y, z])

    if inside_VdW(coords, VdWs, pt):
        return 1

    return 0


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


if __name__ == '__main__':

    args = options()
    struct = np.loadtxt(args.filename, skiprows=2)
    Z_atoms = map(int, struct[:,0].tolist())
    VdWs = [ ELEMENTS[x].vdwrad for x in Z_atoms ]
    coords = struct[:,1:]
    
    vol = calc_vol(coords, VdWs, npts=args.npts)
    
    print(banner(ch="=", length=50))
    print("Molecular Volume (A^3) of %s" % args.filename)
    print("%-18.6f" % vol)
    print
    print("Points used for sampling: %8.2e" % args.npts)
    print(banner(ch="=", length=50))

