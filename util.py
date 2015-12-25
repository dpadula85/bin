#!/usr/bin/env python

import os
import sys
import numpy as np
from itertools import groupby

from elements import ELEMENTS

verbosity = False

# Dictionary for energy conversion
# The key should be of the form u_a2u_b, where u_a and u_b are two different
# units


energy_conversion = {'au'       : 1,
                     'eV'       : 27.21138505,
                     'wn'       : 219474.63,
                     'kj/mol'   : 2625.5,
                     'kcal/mol' : 627.503}

# au to eV = 27.21128505 / 1
# eV to au = 1 / 27.21128505
# u1 to u2 = dict[u2]/dict[u1]


# Dictionaries comparison:
# Options are stored in dictionaries. dictA is the default options dictionary
# and dictB is the one passed to the function.
# We want to compare dictA and dictB. We want to add to dictB all the missing 
# keys in dictA with their value.

def checkfile(filename):

    if not os.path.isfile(filename):
        print(banner(text='ERROR', ch='#', length=80))
        print("File %s not found!" % filename)
        sys.exit()


def dict_compare(dictA, dictB):

    for k in dictA.viewkeys() - dictB.viewkeys():
        dictB[k] = dictA[k]

    return dictB


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    # idxs = idxs.split()

    for idx in idxs:

        to_extend = idx.split('-')

        if len(to_extend) > 1:

            sel =  map(int, to_extend)
            extended += range(sel[0],sel[1]+1,1)

        else:
        
            extended.append(int(idx))
    
    return extended


def compact_extended_list(idxs, factor=0):

    # factor optional parameter to clean out python numeration starting from 0
    compact = []

    for k, iterable in groupby(enumerate(sorted(idxs)), lambda x: x[1] - x[0]):

         rng = list(iterable)

         if len(rng) == 1:

             s = str(rng[0][1] + factor)

         else:

             s = "%s-%s" % (rng[0][1] + factor, rng[-1][1] + factor)

         compact.append(s)

    return compact


def refframe(A, B, C):
    '''Returns a reference frame where the x axis goes from A to B, the y axis
    passes through C and the z axis is built accordingly.'''

    x = (B - A) / np.linalg.norm(B - A)

    # Define the point P on x whose perpendicular to x passes through C
    P = A + np.dot((C - A), x) * x
    y = (C - P) / np.linalg.norm(C - P)

    z = np.cross(x, y)

    ref = np.array([x, y, z])

    return ref


def refframe_var(A, B, C):
    '''Returns a reference frame where A, B and C are contained in the x-y plane,
    A is the origin, the x axis goes from A to B, the y axis is directed towards
    C and the z axis is built accordingly.'''

    x = (B - A) / np.linalg.norm(B - A)

    # Build a temporary y axis to find the real z axis passing through A and
    # perpendicular to the plane defined by A, B and C.
    tmpy = (C - A) / np.linalg.norm(C - A)

    z = np.cross(x, tmpy)
    y = np.cross(z, x)

    ref = np.array([x, y, z])

    return ref


def rot_mat_x(theta):

    theta = np.radians(theta)
    Rx = np.zeros((4,4))
    Rx[0] = np.array([1., 0., 0., 0.])
    Rx[1] = np.array([0., np.cos(theta), -1*np.sin(theta), 0.])
    Rx[2] = np.array([0., np.sin(theta), np.cos(theta), 0.])
    Rx[3] = np.array([0., 0., 0., 1.])

    return Rx


def rot_mat_y(theta):

    theta = np.radians(theta)
    Ry = np.zeros((4,4))
    Ry[0] = np.array([np.cos(theta), 0., np.sin(theta), 0.])
    Ry[1] = np.array([0., 1., 0., 0.])
    Ry[2] = np.array([-1*np.sin(theta), 0., np.cos(theta), 0.])
    Ry[3] = np.array([0., 0., 0., 1.])

    return Ry


def rot_mat_z(theta):

    theta = np.radians(theta)
    Rz = np.zeros((4,4))
    Rz[0] = np.array([np.cos(theta), -1*np.sin(theta), 0., 0.])
    Rz[1] = np.array([np.sin(theta), np.cos(theta), 0., 0.])
    Rz[2] = np.array([0., 0., 1., 0.])
    Rz[3] = np.array([0., 0., 0., 1.])

    return Rz


def rot(axis, theta):
    '''Returns the rotation matrix for the clockwise rotation about
    axis by theta according to Rodrigues' formula.'''

    axis = axis / np.linalg.norm(axis)
    theta = np.radians(theta)

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


def v1v2_angle(v1, v2):

    dotprod = np.dot(v1, v2)
    theta = np.degrees(np.arccos(dotprod / (np.linalg.norm(v1) * np.linalg.norm(v2))))

    return theta


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
                atom[0] = ELEMENTS[atom[0]].symbol
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
    pass
