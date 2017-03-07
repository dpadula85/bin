#!/usr/bin/env python

import os
import sys
import numpy as np
import argparse as arg

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''Calculates the position of a set of points onto
                                            the surface of a sphere.''',
                                            formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Calculation Options
    #
    calc = parser.add_argument_group("Calculation Options")

    calc.add_argument('-r', '--radius', default=1.0, type=float, dest='R',
                     help='''Radius of the Sphere.''')

    calc.add_argument('-m', '--mode', default='rand', type=str, choices=['rand', 'grid', 'fib'],
                      help='''Distribution of the points on the surface of the Sphere.''',
                      dest='Mode')

    calc.add_argument('-n', default=50, type=int, dest='NPts', help='''Number of points
                      randomly distributed on the Sphere.''')

    calc.add_argument('--nphi', default=None, type=int, dest='NPhi',
                      help='''Number of points of the grid along the z axis.''')

    calc.add_argument('--ntheta', default=None, type=int, dest='NTheta',
                      help='''Number of points of the grid along the circle.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--output', default=None, type=str, dest='OutFile',
                     help='''Output File.''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


def sphere_grid(r, nphi, nth):

    lp = nphi * 1j
    lt = nth * 1j
    phi, theta = np.mgrid[0:2*np.pi:lp, 0:np.pi:lt]

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    coords = np.c_[x.ravel(), y.ravel(), z.ravel()]
    coords = unique_rows(coords)

    return coords


def sphere_rand(r, n):

    coords = np.random.normal(size=(n,3))
    coords = unique_rows(coords)
    norm = np.linalg.norm(coords, axis=1)
    coords = np.einsum('ij,i->ij', coords, 1 / norm)

    return coords * r


def sphere_fib(r, n):

    r = float(r)
    angle = np.pi * (3. - np.sqrt(5.))
    theta = np.arange(n) * angle
    step = r / n
    z = np.linspace(r - r / n , r / n - r, n)
    rr = np.sqrt(r**2 - z**2)

    coords = np.zeros((n,3))
    coords[:,0] = rr * np.cos(theta)
    coords[:,1] = rr * np.sin(theta)
    coords[:,2] = z

    return np.array(coords)


def unique_rows(data):

    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))

    return uniq.view(data.dtype).reshape(-1, data.shape[1])


if __name__ == '__main__':

    Opts = options()
    if Opts['Mode'] == 'rand':
        coords = sphere_rand(Opts['R'], Opts['NPts'])

    elif Opts['Mode'] == 'grid':
        coords = sphere_grid(Opts['R'], Opts['NPhi'], Opts['NTheta'])

    elif Opts['Mode'] == 'fib':
        coords = sphere_fib(Opts['R'], Opts['NPts'])

    N = len(coords)
    coords = np.c_[np.ones(N), coords]

    if not Opts['OutFile']:
        Opts['OutFile'] = "sphere_%d" % N

    with open("%s.xyz" % Opts['OutFile'], "w") as f:
        f.write("%d\n\n"% N)
        np.savetxt(f, coords, fmt="%5d %12.6f %12.6f %12.6f")
