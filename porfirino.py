#!/usr/bin/env python

import os
import sys
import numpy as np
import argparse as arg
import shutil as sh

from elements import ELEMENTS

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''
    Projects a reference structure on a set of points along a helix.''')

    # Optional arguments
    parser.add_argument('-n', default=0, type=float, help='''
    Number of centers to project the molecule on.''')

    parser.add_argument('-c', '--centers', default='centers.dat', help='''
    File containing the coordinates of the centers.''')

    parser.add_argument('-r', '--ref', default='ref.inc', help='''
    File containing the reference geometry.''')

    parser.add_argument('-g', '--g09ref', help='''
    G09 optimization for reference geometry.''')

    parser.add_argument('-t', '--theta', default=[0.], nargs='+', type=float, help='''
    Angle for the rotation of the reference structure about its z axis.''')

    parser.add_argument('-a', '--axis', default=['z'], choices=['x', 'y', 'z'], nargs='+',
    type=str, help='''Angle for the rotation of the reference structure about its z axis.''')

    parser.add_argument('-o', '--output', default='final', help='''
    Root name of the .xyz and/or .pdb file to save.''')

    parser.add_argument('--savexyz', action='store_true', default=False, help='''
    Save an .xyz file''')

    parser.add_argument('--savepdb', action='store_true', default=True, help='''
    Save a .pdb file''')

    args = parser.parse_args()

    return args


def checkfile(filename):

    if not os.path.isfile(filename):
        print(banner(text='ERROR', ch='#', length=80))
        print("File %s not found!" % filename)
        sys.exit()


def get_struct(infile):

    structure = []
    opt = 0
    with open(infile, 'r') as f:
        for line in f:
            if "Optimization completed" in line:
                opt = 1
    
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
    
    return np.asarray(structure)


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
            counter_dict = {}

            for atom in molecule:

                i += 1
                atom[0] = ELEMENTS[atom[0]].symbol

                try:
                    counter_dict[atom[0]] += 1

                except KeyError:
                    counter_dict[atom[0]] = 1

                k = counter_dict[atom[0]]
                atom_name = "%s%d" % (atom[0], k)
                f.write(line % (i, atom_name, resname, j, atom[1], atom[2], atom[3], atom[0]))

            # At the end of each molecule
            f.write('TER\n')

        # At the end of the file
        f.write('END')

    return


def write_XYZ(xyzout, coords):

    # Here coords is just an np.array
    line = '%2s %10.6f %10.6f %10.6f'

    with open(xyzout, 'w') as f:

        f.write('%d\n' % len(coords))
        f.write('Title\n')
        np.savetxt(f, coords, fmt=line)

    return


def rot(axis, theta):
    '''Returns the rotation matrix for the anticlockwise rotation about
    axis by theta according to Rodrigues' formula.'''

    axis = axis / np.linalg.norm(axis)
    theta = -1 * np.radians(theta)
    I = np.eye(3)

    # Define axis' cross-product matrix
    K = np.cross(I, axis)

    R = I + np.sin(theta) * K + (1 - np.cos(theta)) * np.linalg.matrix_power(K, 2)

    return R


def transl_mat(v):

    # Define the transformation matrix for a translation
    T = np.eye(4)
    T[-1,:3] = v

    return T


def rototransl(axis, theta, T=np.zeros(3)):
    '''Returns a 4x4 rototranslation matrix, where the rotation part is given
    by the anticlockwise rotation about axis by theta, and the
    translation by the vector T.'''

    R = rot(axis, theta)
    R = np.vstack([R, np.array([0., 0., 0.])])
    R = np.c_[R, np.array([0., 0., 0., 1.])]
    T_mat = transl_mat(T)

    R[-1,:] = T_mat[-1,:]

    return R


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

    print()
    print(banner(ch='=', length=80))
    print(banner(text='porfirino.py', ch=' ', length=80))
    print(banner(ch='=', length=80))
    print()

    args = options()

    cartesian = np.eye(3)
    ux = cartesian[0]
    uy = cartesian[1]
    uz = cartesian[2]

    # Load centers' coordinates
    centers = np.loadtxt(args.centers)

    if args.n == 0:
        n = len(centers)
    
    elif args.n > len(centers):
        print("Not enough centers to process. Only %d available" % len(centers))
        sys.exit()

    else:
        n = args.n


    if args.g09ref:
        # Reference structure from a G09 optimization
        checkfile(args.g09ref)
        opt = get_struct(args.g09ref)

    else:
        # Load reference structure
        checkfile(args.ref)
        opt = np.loadtxt(args.ref)


    # Divide the structure into atoms and coordinates
    struct = opt[:,1:]
    atoms = opt[:,0]
    
    # Stack a column of ones to the structure for
    # dot product with a 4x4 transformation matrix
    struct = np.c_[struct, np.ones(len(struct))]

    # Build the 4x4 rotation matrix for the rotation
    # and rotate the structure
    axis = args.axis
    theta = args.theta
    if len(theta) != len(axis) and len(theta) == 1:
        theta = theta * len(axis)

    for i in range(len(axis)):

        ax = axis[i]
        t = theta[i]

        if ax == 'x':
            ax = ux
        elif ax == 'y':
            ax = uy
        elif ax == 'z':
            ax = uz

        R = rototransl(ax, t)
        struct = np.dot(struct, R)
    
    # Initialize results arrays
    final = np.array([]).reshape(0,4)
    final_structure = []

    # For each center we have to define a local reference frame
    # such that the axes in that point (on the border of the helix)
    # are z', parallel to the helix's axis (z in the global frame),
    # y', radial, and x', tangent.
    # This is because in my reference structure the porphirine plane
    # is in the x-y plane, and when I project it on the border of the helix
    # I want the plane of the porphirine to be orthogonal to the helix's axis,
    # thus the x'-y' plane should be orthogonal to it.
    # Should I need a different arrangement of the porphirine, I may want
    # to change the reference frame on the helix's border (or in the ref structure).
    for center in centers[:n,:]:
    
        # z' is equal to the global z axis. It is only translated
        # by the center's vector
        zc = uz 
        
        # y' is orthogonal to z' and passes through center.
        # To find it, project the center's vector onto z',
        # and define the distance between the center's vector
        # and its projection onto z'.
        A = np.dot(center, zc) * zc
        yc = center - A
        yc = yc / np.linalg.norm(yc) 
    
        # x' is orthogonal to the other two axes, and the reference
        # frame on the border of the helix is completed.
        xc = np.cross(zc, yc)
        ref = np.array([xc, yc, zc])
    
        # Build the 4x4 rotation matrix
        R = np.dot(cartesian, ref.T)
        R = np.vstack([R, np.array([0., 0., 0.])])
        R = np.c_[R, np.array([0., 0., 0., 1.])]
    
        # Build the 4x4 translation matrix
        T = transl_mat(center)
   
        # Transform the coordinates. First rotation, then translation
        transformed = np.dot(struct, R)
        transformed = np.dot(transformed, T)
    
        # Remove the column of ones needed for appropriate dimensions
        # during the transformation and reappend the column of atomic weights
        transformed = np.c_[atoms, transformed[:,:3]]
   
        # Append the structure of this translated molecule to an array for
        # final writing onto a file
        final = np.vstack((final, transformed))
        final_structure.append(transformed.tolist())

    # Here the for cycle is finished!!
    if args.savepdb:
        write_PDB('%s.pdb' % args.output, final_structure)
        print("Output saved in %s.pdb" % args.output)

    if args.savexyz:
        write_XYZ('%s.xyz' % args.output, final)
        print("Output saved in %s.xyz" % args.output)

    #
    # Generation of files for De Voe Calculation
    #
    if os.path.isdir(args.output):
        sh.rmtree(args.output)

    os.makedirs(args.output)

    # Interaction matrix
    with open('matrix.txt', 'w') as f:

        f.write("%d\n" % n)

        for i in range(1, n + 1):
            f.write("1 %5d %5d\n" % (2*i -1, 2*i))

    sh.move('matrix.txt', args.output)

    # Dipoles
    # N atoms are 9, 14, 19, 24
    # meso C atoms are 1, 2, 3, 4
    # The molecule has 92 atoms, 64 without H atoms
    CA, CB, CC, CD = 1, 2, 3, 4 
    natoms = len(struct)
    with open('dipoles.txt', 'w') as f:

        while CA < len(final):

            f.write("%5d %5d 0 0. 0\n" % (CA, CC))
            f.write("%5d %5d 0 0. 1\n" % (CA, CC))
            f.write("%5d %5d 0 0. 1\n" % (CB, CD))

            CA += natoms
            CB += natoms
            CC += natoms
            CD += natoms

    sh.move('dipoles.txt', args.output)

    # Mol file
    np.savetxt('mol.txt', final[:,1:], fmt='%10.6f')
    sh.move('mol.txt', args.output)
    print()
    print(banner(ch='=', length=80))
    print()
