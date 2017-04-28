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

    parser.add_argument('-sf', '--sf', choices=['g09in', 'pdb', 'mol2', 'xyz'], help='''Format of the input structure.''')

    parser.add_argument('-n', '--nunits', type=int, default=1, help='''Number of times to repeat the unit.''')

    parser.add_argument('-d', '--distance', type=float, default=0., help='''Distance between close units.''')

    parser.add_argument('-dd', '--dd', default=False, action='store_true', help='''Increase distance for each unit.''')

    parser.add_argument('-a', '--angle', type=float, default=5., help='''Angle between close units.''')

    parser.add_argument('-ax', '--axis', choices=['x', 'y', 'z'], default='x', help='''Propagation axis for the repetition of the units.''')

    parser.add_argument('-o', '--output', type=str, default='out', help='''Output prefix.''')

    args = parser.parse_args()

    return args


def checkfile(filename):

    if not os.path.isfile(filename):
        print("File %s not found!" % filename)
        sys.exit()


def flatten(lst):
    return sum( ([x] if not isinstance(x, list) else flatten(x) for x in lst), [] )


def parse_G09_input(infile):
    '''Returns the structure.'''

    checkfile(infile)
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


def parse_PDB(pdb_file):

    checkfile(pdb_file)

    atom_names = []
    atoms = []
    res_names = []
    res_ids = []
    atom_coord = []

    with open(pdb_file) as f:

        for line in f:

            # Read Atoms
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':

                atom_names.append(line[12:15])
                atom = line[12:15]
                atom = ''.join([i for i in atom if not i.isdigit()]).strip()
                atoms.append(atom)
                res_names.append(line[17:19])

                try:
                    res_ids.append(int(line[22:25]))

                except ValueError:
                    res_ids.append(None)

                x = float(line[30:37])
                y = float(line[38:45])
                z = float(line[46:53])
                atom_coord.append([x, y, z])

    struct = [ flatten(list(el)) for el in zip(atoms, atom_coord) ]
    return struct


def parse_XYZ(xyz_file):

    checkfile(xyz_file)
    struct = np.genfromtxt(xyz_file, skip_header=2, dtype=None)

    return struct.tolist()


def parse_MOL2(mol2_file):

    checkfile(mol2_file)

    atom_names = []
    atom_types = []
    atoms = []
    res_names = []
    res_ids = []
    atom_coord = []

    with open(mol2_file) as f:

        while True:

            line = f.readline()

            if not line:
                break

            # Read initial lines
            if line[0:17] == '@<TRIPOS>MOLECULE':
                f.readline()
                info = f.readline().split()
                NAtoms = int(info[0])
                try:
                    NRes = int(info[2])
                except:
                    NRes = 1

            # Read Atoms
            elif  line[0:13] == '@<TRIPOS>ATOM':
                for i in range(NAtoms):

                    data = f.readline().split()
                    
                    # save data for the old one
                    atom_names.append(data[1])
                    atom = ''.join([i for i in data[1] if not i.isdigit()]).strip()
                    atoms.append(atom)
                    atom_types.append(data[5])
                    atom_coord.append(map(float, data[2:5]))

    struct = [ flatten(list(el)) for el in zip(atoms, atom_coord) ]
    return struct


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
    sf = args.sf
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

    #
    # Get structure
    #
    if sf == 'mol2':
        m = parse_MOL2(infile)

    elif sf == 'pdb':
        m = parse_PDB(infile)

    elif sf == 'xyz':
        m = parse_XYZ(infile)

    elif sf == 'g09in':
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
