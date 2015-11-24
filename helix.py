#!/usr/bin/env python

import numpy as np
import argparse as arg

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


def rot_mat_z(theta):

    theta = np.radians(theta)
    Rz = np.zeros((4,4))
    Rz[0] = np.array([np.cos(theta), -1*np.sin(theta), 0., 0.])
    Rz[1] = np.array([np.sin(theta), np.cos(theta), 0., 0.])
    Rz[2] = np.array([0., 0., 1., 0.])
    Rz[3] = np.array([0., 0., 0., 1.])

    return Rz


def v1v2_angle(v1, v2):

    dotprod = np.dot(v1, v2)
    theta = np.degrees(np.arccos(dotprod / (np.linalg.norm(v1) * np.linalg.norm(v2))))

    return theta


def dist(p1, p2, normal=None):
    '''Returns a tuple containing the distance between two points
    and the distance on the plane whose normal is supplied.'''

    # Calculate projection on to xy plane by default
    if not normal:
        normal = np.array([0., 0., 1.])

    d = p2 - p1
    dplane = d - np.dot(d, normal)*normal
    theta = v1v2_angle(d, dplane)

    return np.linalg.norm(d), np.linalg.norm(dplane), theta


def write_PDB(pdbout, coords):

    # For better organization of the output writing
    # coords must be a list of lists:
    # coords = [[at1mol1, at2mol1, ...], [at1mol2, at2mol2, ...], ..., [at1molN, at2molN, ...]]

    fmt = "ATOM  %5d %-4s %3s %5d    %8.3f%8.3f%8.3f  0.00  0.00  %s\n"
    resname = 'PRF'

    with open(pdbout, 'w') as f:

        # i : total atom counter
        # j : residue counter
        # k : atom in molecule counter
        i = 0
        j = 0

        for molecule in coords:

            j += 1
            k = 0

            for atom in molecule:

                k += 1
                i += 1
                atom[0] = ELEMENTS[atom[0]].symbol
                atom_name = "%s%d" % (atom[0], k)
                f.write(fmt % (i, atom_name, resname, j, atom[1], atom[2], atom[3], atom[0]))

    return


def write_XYZ(xyzout, coords):

    # Here coords is just an np.array
    line = '%2s %10.6f %10.6f %10.6f'

    with open(xyzout, 'w') as f:

        f.write('%d\n' % len(coords))
        f.write('Title\n')
        np.savetxt(f, coords, fmt=line)

    return


if __name__ == '__main__':

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

    if hs == 'r':
        factor = -1
    elif hs == 'l':
        factor = 1

    cartesian = np.eye(3)
    ux = cartesian[0]
    uy = cartesian[1]
    uz = cartesian[2]

    # Generate S helix and save its structure
    s_helix = helix_points(r, 1, ss, ps, hs)
    s_helix = s_helix[:-2]
    # s_helix = ellipse_points(r, turns, ss, ps, 5, 'l')
    # s_helix = circle_points(r, ps)
    s_helix_save = np.c_[np.ones(len(s_helix)), s_helix]
    write_XYZ('helix_%d_%s.xyz' % (ps, hs), s_helix_save)

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
        Rz = rot_mat_z(angle)
        T = np.eye(4)
        T[-1,:3] = np.array([0., 0., point[2]])

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
    write_XYZ('helix_%d_%s.xyz' % (pb, hb), b_helix_save)

    # final = final[final[:,2].argsort()]
    final_save = np.c_[np.ones(len(final)), final]
    write_XYZ('%s.xyz' % args.output, final_save)
    print("Output saved in %s" % args.output)
