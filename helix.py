#!/usr/bin/env python

import numpy as np
import argparse as arg

import util as u

# Parameters that define a helix
# r : radius of the helix
# turns : number of turns to compute
# ppt : points per turn
# s : step of the helix

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''
    Calculates the coordinates of points on a helix.''')

    # Optional arguments
    parser.add_argument('-r', '--radius', default=80, type=float, help='''
    Helix Radius.''')

    parser.add_argument('-t', '--turns', default=2, type=float, help='''
    Number of turns of the helix.''')

    parser.add_argument('--ss', default=55, type=float, help='''
    Step of S helix (distance between two turns).''')

    parser.add_argument('--sb', default=165, type=float, help='''
    Step of B helix (distance between two turns).''')

    parser.add_argument('--ratio', default=3, type=float, help='''
    Ratio between the steps of B and S helices.''')

    parser.add_argument('--ps', default=24, type=float, help='''
    Points per turn to calculate for S helix.''')

    parser.add_argument('--pb', default=72, type=float, help='''
    Points per turn to calculate for B helix.''')

    parser.add_argument('--skew', default=0, type=float, help='''
    Skew angle for S helix.''')

    parser.add_argument('--flat', action='store_true', default=False,
    help='''Flatten S helix to a circle.''')

    parser.add_argument('-o', '--output', default='tube', help='''
    Root of the name of the output file.''')

    parser.add_argument('--hs', choices=['r', 'l'], default='r', help='''
    Helicity of S helix.''')

    parser.add_argument('--hb', choices=['r', 'l'], default='r', help='''
    Helicity of B helix.''')

    args = parser.parse_args()

    return args


def helix_points(radius, turns, step, ppt, helicity='r'):

    angles = np.arange(0, 360 * turns + 1, 360 / ppt) 
    points = np.array([]).reshape(0,3)
    # points = np.vstack((points, np.array([0., 0., 0.])))

    if helicity == 'r':
        helicity = 1

    elif helicity == 'l':
        helicity = -1

    for t in angles:

        x = radius * np.cos(np.radians(t))
        y = radius * np.sin(np.radians(t)) * helicity
        z = step * np.radians(t) / (2*np.pi)
        p = np.array([x, y, z])
        points = np.vstack((points, p))

    return points


def ellipse_points(radius, turns, step, ppt, skew_angle=0., helicity='r'):

    # When skew_angle = 0, ellipse_points = helix_points
    angles = np.arange(0, 360 * turns, 360 / ppt) 
    points = np.array([]).reshape(0,3)
    # points = np.vstack((points, np.array([0., 0., 0.])))

    if helicity == 'r':
        helicity = 1

    elif helicity == 'l':
        helicity = -1

    for t in angles:

        x = radius * np.cos(np.radians(t))
        y = radius * np.sin(np.radians(t)) * helicity
        z = step * np.radians(t) / (2*np.pi) + radius * np.sin(np.radians(t)) * np.tan(np.radians(skew_angle)) 
        p = np.array([x, y, z])
        points = np.vstack((points, p))

    return points


def circle_points(radius, ppt):

    angles = np.arange(0, 360, 360 / ppt) 
    points = np.array([]).reshape(0,3)
    # points = np.vstack((points, np.array([0., 0., 0.])))

    for t in angles:

        x = radius * np.cos(np.radians(t))
        y = radius * np.sin(np.radians(t))
        z = 0
        p = np.array([x, y, z])
        points = np.vstack((points, p))

    return points


def dist(p1, p2, normal=None):
    '''Returns a tuple containing the distance between two points
    and the distance on the plane whose normal is supplied.'''

    # Calculate projection on to xy plane by default
    if not normal:
        normal = np.array([0., 0., 1.])

    d = p2 - p1
    dplane = d - np.dot(d, normal)*normal
    theta = u.v1v2_angle(d, dplane)

    return np.linalg.norm(d), np.linalg.norm(dplane), theta


if __name__ == '__main__':

    print
    print u.banner(ch='=', length=80)
    print u.banner(text='helix.py', ch=' ', length=80)
    print u.banner(ch='=', length=80)
    print

    args = options()

    r = args.radius
    turns = args.turns
    ss = args.ss
    sb = args.sb
    ps = args.ps
    pb = args.pb
    ratio = args.ratio
    hs = args.hs
    hb = args.hb
    skew = args.skew

    if hs == 'r':
        factor = -1
    elif hs == 'l':
        factor = 1

    cartesian = np.eye(3)
    ux = cartesian[0]
    uy = cartesian[1]
    uz = cartesian[2]

    # Generate S helix and save its structure

    if args.flat:
        s_helix = circle_points(r, ps)

    else:
        s_helix = ellipse_points(r, 1.05, ss, ps, skew, hs)

    if hs == 'r':
        s_helix = s_helix[:ps + 2]

    elif hs == 'l':
        s_helix = s_helix[:ps -2]

    # s_helix = ellipse_points(r, turns, ss, ps, 5, 'l')
    # s_helix = circle_points(r, ps)
    s_helix_save = np.c_[np.ones(len(s_helix)), s_helix]

    u.write_XYZ('helix_S%d%s.xyz' % (ps, hs), s_helix_save)
    print("S helix saved in helix_S%d%s.xyz" % (ps, hs))

    # Take two adjacent points on this helix and calculate
    # their distance and their distance onto the x-y plane.
    ps1 = s_helix[0]
    ps2 = s_helix[1]
    ds, dsxy, theta = dist(ps1, ps2)
    print("S Helix: d = %5.2f A  dxy = %5.2f A  theta = %6.2f" % (ds, dsxy, theta))

    # Initialize the array to store B helices
    final = np.array([]).reshape(0,3)

    i = 0
    for point in s_helix:

        # I gotta generate a helix starting in point
        # To do this, I generate the helix in the global reference
        # frame, and then I rotate it about z of 2*pi/ps and translate
        # along z by point's z component
        angle = np.degrees(factor*i*2*np.pi/ps)
        b_helix = helix_points(r, turns, sb, pb, hb)

        # Stack a column of ones for dot product with 4x4 transformation
        # matrix
        b_helix = np.c_[b_helix, np.ones(len(b_helix))]

        # Generate 4x4 transformation matrices
        Rz = u.rot_mat_z(angle)
        T = u.transl_mat(np.array([0., 0., point[2]]))

        # Transform the helix
        transformed = np.dot(b_helix, Rz)
        transformed = np.dot(transformed, T)

        # Remove the column of ones needed for appropriate dimensions
        # during the transformation and reappend the column of atomic weights
        transformed = transformed[:,:3]

        final = np.vstack((final, transformed))
        i += 1
   
    # Take two adjacent points on B helix and calculate
    # their distance and their distance onto the x-y plane.
    pb1 = final[0]
    pb2 = final[1]
    db, dbxy, theta = dist(pb1, pb2)
    print("B Helix: d = %5.2f A  dxy = %5.2f A  theta = %6.2f" % (db, dbxy, theta))

    # Sort points coordinates by z component value
    # for easier selection of a small portion of the tube
    # in VMD
    b_helix_save = np.c_[np.ones(len(final[:pb*turns])), final[:pb*turns]]
    u.write_XYZ('helix_B%d%s.xyz' % (pb, hb), b_helix_save)
    print("B helix saved in helix_B%d%s.xyz" % (pb, hb))
    print

    # final = final[final[:,2].argsort()]
    final_save = np.c_[np.ones(len(final)), final]
    u.write_XYZ('%s_S%d%s_B%d%s.xyz' % (args.output, ps, hs, pb, hb), final_save)
    print("Output saved in %s_S%d%s_B%d%s" % (args.output, ps, hs, pb, hb))

    print
    print(u.banner(ch='=', length=80))
    print
