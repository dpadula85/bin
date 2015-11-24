#!/usr/bin/env python

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

def dict_compare(dictA, dictB):

    for k in dictA.viewkeys() - dictB.viewkeys():
        dictB[k] = dictA[k]

    return dictB


def rot_mat_x(theta):

    theta = np.radians(theta)
    Rx = np.zeros((4,4))
    Rx[0] = np.array([1., 0., 0., 0.])
    Rx[1] = np.array([0., np.cos(theta), -1*np.sin(theta), 0.])
    Rx[2] = np.array([0.,cnp.sin(theta), np.cos(theta), 0.])
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


def v1v2_angle(v1, v2):

    dotprod = np.dot(v1, v2)
    theta = np.degrees(np.arccos(dotprod / (np.linalg.norm(v1) * np.linalg.norm(v2))))

    return theta


def write_PDB(pdbout, coords):

    # For better organization of the output writing
    # coords must be a list of lists:
    # coords = [[at1mol1, at2mol1, ...], [at1mol2, at2mol2, ...], ..., [at1molN, at2molN, ...]]

    fmt = "ATOM  %5d %-4s %3s %5d    %8.3f%8.3f%8.3f  0.00  0.00  %s\n"
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
                f.write(fmt % (i, atom_name, resname, j, atom[1], atom[2], atom[3], atom[0]))

            # At the end of each molecule
            f.write('TER')

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
