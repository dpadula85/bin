#!/usr/bin/env python

import os
import sys
import numpy as np
from itertools import groupby

from elements import ELEMENTS

verbosity = False

# Dictionary for energy conversion


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


def flatten(lst):
    return sum( ([x] if not isinstance(x, list) else flatten(x) for x in lst), [] )


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


def symm_mat(M):
    '''Symmetrize an upper- or lower diagonal matrix.'''
    return M + M.T - np.diag(M.diagonal())


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


def rototransl(axis, theta, T=None)):
    '''Returns a 4x4 rototranslation matrix, where the rotation part is given
    by the anticlockwise rotation about axis by theta, and the
    translation by the vector T.'''

    R = rot(axis, theta)
    R = np.vstack([R, np.array([0., 0., 0.])])
    R = np.c_[R, np.array([0., 0., 0., 1.])]
    if T:
        T_mat = transl_mat(T)

    R[-1,:] = T_mat[-1,:]

    return R


def v1v2_angle(v1, v2):
    '''Returns the angle between two vectors.'''

    dotprod = np.dot(v1, v2)
    theta = np.degrees(np.arccos(dotprod / (np.linalg.norm(v1) * np.linalg.norm(v2))))

    return theta


def angle(A, B, C):
    '''Returns the bond angle defined by atoms ABC.'''

    ab = B - A
    bc = C - B

    return v1v2_angle(ab, bc)


def dihedral(A, B, C, D):
    '''Returns the dihedral angle between the planes containing bonds AB and CD.'''

    ab = B - A
    bc = C - B
    cd = D - C
    
    c1 = np.cross(ab, bc)
    n1 = c1 / np.linalg.norm(c1)
    c2 = np.cross(bc, cd)
    n2 = c2 / np.linalg.norm(c2)
    m1 = np.cross(n1, bc)
    x = np.dot(n1, n2) 
    y = np.dot(m1, n2)

    return np.arctan2(y, x) * 180 / np.pi


def write_PDB(pdbout, coords):

    # For better organization of the output writing
    # coords must be a list of lists of lists:
    # coords = [[at1mol1, at2mol1, ...], [at1mol2, at2mol2, ...], ..., [at1molN, at2molN, ...]],
    # where atXmolY is a list of four elements: symbol (or atomic number) and x, y, z coordinates

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


def parse_MOL2(mol2file):
    '''The following code creates this structure:
    res_names and res_ids are plain lists, and the elements
    appear only once per residue.
    atom_names, type and coord are lists of lists. Each sublist
    corresponds to a residue and contains all the info
    for  on the atoms in that residue.'''


    with open(mol2file) as f:

        FoundAt = False
        FoundBond = False

        while True:
            line = f.readline()
            if not line:
                break

            # skip comments
            elif line[0] == '#':
                continue
            
            # Read initial lines
            elif line[0:17] == '@<TRIPOS>MOLECULE':
                f.readline()
                info = f.readline().split()
                NAtoms = int(info[0])
                NBonds = int(info[1])
                try:
                    NRes = int(info[2])
                except:
                    NRes = 1

                atom_names = []
                atom_types = []
                res_names = []
                res_ids = []
                atom_coord = []
                Ib1 = []
                Ib2 = []

            # Read Atoms
            elif  line[0:13] == '@<TRIPOS>ATOM':
                for i in range(NAtoms):
                    data = f.readline().split()
                    
                    # Special case for the first residue
                    try:
                        res_ids[-1]
                    except IndexError:
                        res_ids.append(int(data[6]))
                    
                    # Special case for the first residue
                    try:
                        res_names[-1]
                    except IndexError:
                        res_names.append(data[7])
                    
                    # Check id: new residue or old one
                    # if in new residue
                    if int(data[6]) != res_ids[-1]:
                        res_ids.append(int(data[6]))
                        res_names.append(data[7])
                        
                        # save data for the old one
                        atom_names.append(names)
                        atom_types.append(types)
                        atom_coord.append(coords)
                    
                        # initialize data for the new one
                        names = [data[1]]
                        types = [data[5]]
                        # coords = [data[2:5]]
                        coords = map(float, data[2:5])
                        coords = [coords]
                    
                    # if still in the old residue
                    else:
                    
                        # save data for the new atom
                        try:
                            names.append(data[1])
                            types.append(data[5])
                            coords.append(map(float, data[2:5]))
                    
                        # unless no atom has been saved before
                        except:
                            names = [data[1]]
                            types = [data[5]]
                            # coords = [data[2:5]]
                            coords = map(float, data[2:5])
                            coords = [coords]

                FoundAt = True
                if FoundAt:
                    # save data for the last residue 
                    atom_names.append(names)
                    atom_types.append(types)
                    atom_coord.append(coords)

            elif line[0:13] == '@<TRIPOS>BOND':
                for i in range(NBonds):
                    data = f.readline().split()[1:3]
                    Ib1 += [int(data[0])]
                    Ib2 += [int(data[1])]

                FoundBond = True

            elif  line[0:21] == '@<TRIPOS>SUBSTRUCTURE':
                if FoundAt and FoundBond:
                    break 

    # Build Connectivity Matrix
    Conn = np.zeros((NAtoms, 9), dtype=int)
    for i in range(NBonds):
        Idx1 = np.argmax(Conn[Ib1[i]-1] == 0)
        Idx2 = np.argmax(Conn[Ib2[i]-1] == 0)
        Conn[Ib1[i]-1][Idx1] = Ib2[i]
        Conn[Ib2[i]-1][Idx2] = Ib1[i]

    return atom_names, atom_types, res_names, res_ids, atom_coord, Conn


def get_group(D, connectivity, visited=None, C=None):
    '''Returns the list of the indexes of the group of atoms connected
    to atom D, given the connectivity matrix. Atom C is optional and
    should be linked to D as well. C is useful to get all the atoms
    connected to D excluding the portion in which C is included.
    However, atom C is included in the list returned.
    This function is useful to define the groups of atoms which should
    undergo a certain transformation, for example for a dihedral scan.'''

    if not visited:
        if C:
            visited = [C, D]
        else:
            visited = [D]

    for atom in connectivity[D]:

        if atom != 0 and atom - 1 not in visited:
            visited.append(atom - 1)
            get_group(atom - 1, connectivity, visited)
            
        if atom - 1 in visited:
            continue

    return sorted(visited)


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
