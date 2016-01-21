#!/usr/bin/env python

import re
import os
import sys
import numpy as np
import argparse as arg


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Repeats a monomeric unit in space.')

    # Optional arguments
    parser.add_argument('-f', '--filename', help='''Input structure.''')

    parser.add_argument('-sf', '--format', help='''Format of the input structure.''')

    parser.add_argument('-n', '--nunits', type=int, default=1, help='''Number of times to repeat the unit.''')

    parser.add_argument('-d', '--distance', type=float, default=4., help='''Distance between close units.''')

    parser.add_argument('-dd', '--dd', default=False, action='store_true', help='''Increase distance for each unit.''')

    parser.add_argument('-a', '--angle', type=float, default=5., help='''Angle between close units.''')

    parser.add_argument('-ax', '--axis', choices=['x', 'y', 'z'], help='''Propagation axis for the repetition of the units.''')

    parser.add_argument('-o', '--output', type=str, default='out', help='''Output prefix.''')

    args = parser.parse_args()

    return args


def parse_G09_input(infile):
    '''Returns the structure.'''

    structure = []
    opt = 0

    # Pattern describing the charge-multiplicity line
    pattern = re.compile('(\d\s\d)')

    with open(infile, 'r') as f:
        for line in f:
            if pattern.match(line):
                opt = 1

            if not line.strip():
                opt = 0

            if opt == 1 and not pattern.match(line):
                curr_line = filter(None, line.split())
                atom = curr_line[0]
                atom_x = float(curr_line[1])
                atom_y = float(curr_line[2])
                atom_z = float(curr_line[3])
                structure.append([atom, atom_x, atom_y, atom_z])

    return structure


def write_PDB(pdbout, coords):

    # For better organization of the output writing
    # coords must be a list of lists:
    # coords = [[at1mol1, at2mol1, ...], [at1mol2, at2mol2, ...], ..., [at1molN, at2molN, ...]]

    line = "ATOM  %5d %-4s %3s %5d    %8.3f%8.3f%8.3f  0.00  0.00  %s\n"
    resname = 'MOL'

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
                # atom[0] = ELEMENTS[atom[0]].symbol
                atom_name = "%s%d" % (atom[0], k)
                f.write(line % (i, atom_name, resname, j, atom[1], atom[2], atom[3], atom[0]))

            # At the end of each molecule
            f.write('TER\n')

        # At the end of the file
        f.write('END')

    return


def rot(axis, theta):
    '''Returns the rotation matrix for the anticlockwise rotation about
    axis by theta according to Rodrigues' formula.'''

    axis = axis / np.linalg.norm(axis)
    theta = -1 * np.radians(theta)

    # Define axis' cross-product matrix
    K = np.zeros((3,3))
    K[0] = [       0, -axis[2],  axis[1]]
    K[1] = [ axis[2],        0, -axis[0]]
    K[2] = [-axis[1],  axis[0],        0]

    I = np.eye(3)
    R = I + np.sin(theta) * K + (1 - np.cos(theta)) * np.linalg.matrix_power(K, 2)

    return R


def transl_mat(v):

    # Define the transformation matrix for a translation
    T = np.eye(4)
    T[-1,:3] = v

    return T


def rototransl(axis, theta, T):
    '''Returns a 4x4 rototranslation matrix, where the rotation part is given
    by the clockwise rotation about axis by theta, and the
    translation by the vector T.'''

    R = rot(axis, theta)
    R = np.vstack([R, np.array([0., 0., 0.])])
    R = np.c_[R, np.array([0., 0., 0., 1.])]
    T_mat = transl_mat(T)

    R[-1,:] = T_mat[-1,:]

    return R


if __name__ == '__main__':

    args = options()
    infile = args.filename
    nunits = args.nunits
    distance = args.distance
    angle = args.angle
    axis = args.axis
    outfile = args.output

    cartesian = np.eye(3)

    if axis == 'x':
        axis = cartesian[0]

    elif axis == 'y':
        axis = cartesian[1]

    elif axis == 'z':
        axis = cartesian[2]

    m = parse_G09_input(infile)
    natoms = len(m)
    struct = np.array([ [atom[1], atom[2], atom[3]] for atom in m ])
    atoms = [ atom[0] for atom in m ]
    struct1 = np.c_[struct, np.ones(len(struct))]

    for idx,i in enumerate(range(nunits), start=1):

        if args.dd:
            d = idx * distance
        else:
            d = distance

        a = idx * angle

        M = rototransl(axis, a, d * axis)
        transformed = np.dot(struct1, M)
        transformed = transformed[:,:3]
        
        struct = np.vstack((struct, transformed))

        step_struct = [[ [ el[0], el[1][0], el[1][1], el[1][2] ] for el in zip(atoms, transformed.tolist()) ]]
        step_out = '%s_%dA_%ddeg.pdb' % (outfile, d, a)
        write_PDB(step_out, step_struct)


    atoms = np.tile(atoms, nunits + 1).tolist()
    final_struct = struct.tolist()

    data = [ [ el[0], el[1][0], el[1][1], el[1][2] ] for el in zip(atoms, final_struct) ]
    final = [ data[x:x + natoms] for x in range(0, len(data), natoms) ]

    write_PDB('%s_total.pdb' % outfile, final)
