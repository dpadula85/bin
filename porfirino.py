#!/usr/bin/env python

import os
import sys
import numpy as np
import argparse as arg
import shutil as sh

from elements import ELEMENTS
import util as u

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

    parser.add_argument('-t', '--theta', default=0, type=float, help='''
    Angle for the rotation of the reference structure about its z axis.''')

    parser.add_argument('-o', '--output', default='final', help='''
    Root name of the .xyz and/or .pdb file to save.''')

    parser.add_argument('--savexyz', action='store_true', default=False, help='''
    Save an .xyz file''')

    parser.add_argument('--savepdb', action='store_true', default=True, help='''
    Save a .pdb file''')

    args = parser.parse_args()

    return args


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


if __name__ == '__main__':

    print
    print u.banner(ch='=', length=80)
    print u.banner(text=sys.argv[0], ch='=', length=80)
    print u.banner(ch='=', length=80)
    print

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
        u.checkfile(args.g09ref)
        opt = get_struct(args.g09ref)

    else:
        # Load reference structure
        u.checkfile(args.ref)
        opt = np.loadtxt(args.ref)


    # Divide the structure into atoms and coordinates
    struct = opt[:,1:]
    atoms = opt[:,0]
    
    # Stack a column of ones to the structure for
    # dot product with a 4x4 transformation matrix
    struct = np.c_[struct, np.ones(len(struct))]

    # Build the 4x4 rotation matrix for the rotation
    # of the reference structure about the z-axis
    theta = args.theta
    R = u.rot_mat_y(theta) 
    # R = u.rot_mat_z(theta) 

    # Rotation of the reference structure about its z axis.
    # Default is 0 degrees.
    struct = np.dot(struct, R)
    
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
        T = u.transl_mat(center)
   
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
        u.write_PDB('%s.pdb' % args.output, final_structure)
        print("Output saved in %s.pdb" % args.output)

    if args.savexyz:
        u.write_XYZ('%s.xyz' % args.output, final)
        print("Output saved in %s.xyz" % args.output)

    #
    # Generation of files for De Voe Calculation
    #
    if os.path.isdir('DeVoe'):
        sh.rmtree('DeVoe')

    os.makedirs('DeVoe')

    # Interaction matrix
    with open('matrix.txt', 'w') as f:

        f.write("%d\n" % n)

        for i in range(1, n + 1):
            f.write("1 %5d %5d\n" % (2*i -1, 2*i))

    sh.move('matrix.txt', 'DeVoe')

    # Dipoles
    # N atoms are 9, 14, 19, 24
    # The molecule has 92 atoms
    NA, NB, NC, ND = 9, 14, 19, 24
    natoms = len(struct)
    with open('dipoles.txt', 'w') as f:

        while NA < len(final):

            f.write("%5d %5d 0 0. 0\n" % (NA, NC))
            f.write("%5d %5d 0 0. 1\n" % (NA, NC))
            f.write("%5d %5d 0 0. 1\n" % (NB, ND))

            NA += natoms
            NB += natoms
            NC += natoms
            ND += natoms

    sh.move('dipoles.txt', 'DeVoe')

    # Mol file
    np.savetxt('mol.txt', final[:,1:], fmt='%10.6f')
    sh.move('mol.txt', 'DeVoe')

    print
    print u.banner(ch='=', length=80)
    print