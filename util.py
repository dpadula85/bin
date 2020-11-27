#!/usr/bin/env python

import numpy as np
import multiprocessing as mp
from itertools import groupby
from scipy.spatial.distance import cdist


def skiplines(openfile, nlines=0):
    '''
    Function to skip nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.

    Parameters
    ----------
    openfile: object.
        File object to process.
    nlines: int.
        Number of lines to skip.

    Returns
    -------
    line: string.
        Line after skipping nlines + 1 lines.
    '''

    for i in range(nlines + 1):
        line = next(openfile)

    return line


def parallel_fn(fn, iterable, nproc=None):

    '''
    Function to execute a generic function in parallel applying it to all
    the elements of an iterable. The function fn should contain appropriate
    error handling to avoid mp.Pool to hang up.

    Parameters
    ----------
    fn: object.
        Python standalone function.
    iterable: iterable.
        Collection of elements on which fn should be applied.
    nproc: integer (default: None).
        Number of processors for the parallelisation. If not specified all
        available processors will be used.

    Returns
    -------
    data: list.
        List of the returns of function fn for each element of iterable.
    '''

    if not nproc:
        nproc = os.cpu_count()

    pool = mp.Pool(nproc)
    data = pool.map(fn, iterable)                                                                                                                     
    pool.close()
    pool.join()

    return data


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    idxs = idxs.split()

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
    '''
    Recursive function to flatten a nested list.

    Parameters
    ----------
    lst: list.
        Nested list to be flattened.

    Returns
    -------
    flattened: list.
        Flattened list.
    '''

    flattened = sum( ([x] if not isinstance(x, list)
                     else flatten(x) for x in lst), [] )

    return flattened


def centroid(coords, masses=None):
    '''
    Function to compute the centre (or the centre of mass) of a set of
    coordinates.

    Parameters
    ----------
    coord: np.array (N,3).
        coordinates.
    masses: np.array (N) (default: None).
        masses.

    Returns
    -------
    com: np.array (3).
        centre (or centre of mass) of the set of coordinates.
    '''

    com = np.average(coords, axis=0, weights=masses)

    return com


def symm_mat(M):
    '''
    Function to symmetrize an upper- or lower diagonal matrix.

    Parameters
    ----------
    M: np.array (N,N).
        Matrix to be symmetrised.

    Returns
    -------
    M: np.array (N,N).
        Symmetrised matrix.
    '''

    M = M + M.T - np.diag(M.diagonal())

    return M


def get_connectivity_matrix(coords, radii):
    '''
    Function to symmetrize an upper- or lower diagonal matrix.

    Parameters
    ----------
    coords: np.array (N,3).
        Coordinates matrix.
    radii: np.array (N).
        Radii vector.

    Returns
    -------
    M: np.array (N,P).
        Connectivity matrix.
    '''

    # Compute distance matrix
    D = cdist(coords, coords)

    # Add fictitious distance on the diagonal to not find the atom itself
    D = D + np.diag(D.diagonal() + 10)

    # Compute radii sum matrix
    R = radii.reshape(1, -1) / 2 + radii.reshape(-1, 1) / 2

    # Compare them to figure out connectivity
    idxs = D <= R

    # For each line get True idxs
    connidxs = [ np.where(i == True)[0].tolist() for i in idxs ]
    pad = len(max(connidxs, key=len))
    M = np.array([i + [-1]*(pad - len(i)) for i in connidxs])

    return M


def rot(axis, theta):
    '''
    Returns the rotation matrix for the anticlockwise rotation about an
    arbitrary axis by theta according to Rodrigues' formula.

    Parameters
    ----------
    axis: np.array (3).
        Unit vector describing the rotation axis.
    theta: float.
        Angle of rotation (in degrees).

    Returns
    -------
    R: np.array (3,3).
        Rotation matrix.
    '''

    axis = axis / np.linalg.norm(axis)
    theta = -np.radians(theta)
    I = np.eye(3)

    # Define axis' cross-product matrix
    K = np.cross(I, axis)

    R = I + np.sin(theta) * K + (1 - np.cos(theta)) * np.linalg.matrix_power(K, 2)

    return R


def v1v2_angle(v1, v2):
    '''
    Function to compute the angle between two vectors.

    Parameters
    ----------
    v1: np.array (3).
        First vector.
    v2: np.array (3).
        Second vector.

    Returns
    -------
    theta: float.
        Angle between the vector (in degrees).
    '''

    try:
        theta = np.degrees(np.arccos(
                    np.dot(v1, v2) / ( np.linalg.norm(v1) * np.linalg.norm(v2) )
                ))
    except:
        theta = 0.0

    return theta


def angle(A, B, C):
    '''
    Function to compute the angle defined by points A, B, and C.

    Parameters
    ----------
    A: np.array (2 or 3).
        First point.
    B: np.array (2 or 3).
        Second point, vertex of the angle.
    C: np.array (2 or 3).
        Third point.

    Returns
    -------
    abc: float.
        Angle ABC (in degrees).
    '''

    ab = B - A
    bc = C - B
    abc = v1v2_angle(ab, bc)

    return abc


def dihedral(A, B, C, D):
    '''
    Function to compute the dihedral angle between the planes containing
    segments AB and CD.

    Parameters
    ----------
    A: np.array (3).
        First point.
    B: np.array (3).
        Second point.
    C: np.array (3).
        Third point.
    D: np.array (3).
        Fourth point.

    Returns
    -------
    dihedral: float.
        Dihedral angle defined by AB and CD (in degrees).
    '''

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
    dihedral = np.arctan2(y, x) * 180 / np.pi

    return dihedral


def acf(series):
    '''
    Function to compute the autocorrelation function of a time series.

    Parameters
    ----------
    series: np.array.
        Time series.

    Returns
    -------
    acf_t: np.array.
        Autocorrelation function of the time series.
    '''

    N = len(series)
    avg = np.mean(series)
    c0 = np.sum((series - avg)**2) / N

    def r(j):
        return np.sum((series[:N - j] - avg) * (series[j:] - avg)) / (N - j)

    t = np.arange(N)
    acf_t = map(r, t) / c0

    return acf_t


def parse_MOL2(mol2file, cdim=4):
    '''
    Function to parse a Sybyl MOL2 file.

    Parameters
    ----------
    mol2file: string.
        File to parse
    cdim: integer (default: 4)
        Second dimension of the connectivity matrix, representing the maximum
        number of atoms that can be bonded to a certain atom.

    Returns
    -------
    atom_names: list of lists.
        List of atom names. Each residue is contained in a sublist.
    atom_types: list of lists.
        List of atom typess. Each residue is contained in a sublist.
    res_names:
        List of residue names.
    res_ids:
        List of residue ids.
    atom_coord:
        Array of atom coordinates. Each residue is contained in a subarray.
    conn: np.array (N,cdim).
        connectivity matrix of the system.
    '''

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
                        atom_coord.append(np.array(coords))
                    
                        # initialize data for the new one
                        names = [ data[1] ]
                        types = [ data[5] ]
                        coords = list(map(float, data[2:5]))
                        coords = [ coords ]
                    
                    # if still in the old residue
                    else:
                    
                        # save data for the new atom
                        try:
                            names.append(data[1])
                            types.append(data[5])
                            coords.append(list(map(float, data[2:5])))
                    
                        # unless no atom has been saved before
                        except:
                            names = [ data[1] ]
                            types = [ data[5] ]
                            coords = list(map(float, data[2:5]))
                            coords = [ coords ]

                FoundAt = True
                if FoundAt:
                    # save data for the last residue 
                    atom_names.append(names)
                    atom_types.append(types)
                    atom_coord.append(np.array(coords))

            elif line[0:13] == '@<TRIPOS>BOND':
                for i in range(NBonds):
                    data = f.readline().split()[1:3]
                    Ib1 += [ int(data[0]) ]
                    Ib2 += [ int(data[1]) ]

                FoundBond = True

            elif  line[0:21] == '@<TRIPOS>SUBSTRUCTURE':
                if FoundAt and FoundBond:
                    break 

    atom_coord = np.array(atom_coord)

    # Build Connectivity Matrix
    conn = np.zeros((NAtoms, cdim), dtype=int)
    for i in range(NBonds):
        Idx1 = np.argmax(conn[Ib1[i] - 1] == 0)
        Idx2 = np.argmax(conn[Ib2[i] - 1] == 0)
        conn[Ib1[i] - 1,Idx1] = Ib2[i]
        conn[Ib2[i] - 1,Idx2] = Ib1[i]

    return atom_names, atom_types, res_names, res_ids, atom_coord, conn


def parse_PDB(pdbfile, cdim=4):
    '''
    Function to parse a Protein Data Bank file.

    Parameters
    ----------
    pdbfile: string.
        File to parse
    cdim: integer (default: 4)
        Second dimension of the connectivity matrix, representing the maximum
        number of atoms that can be bonded to a certain atom.

    Returns
    -------
    atom_idxs: list of lists.
        List of atom indexes. Each residue is contained in a sublist.
    atom_names: list of lists.
        List of atom names. Each residue is contained in a sublist.
    atom_types: list of lists.
        List of atom typess. Each residue is contained in a sublist.
    res_names:
        List of residue names.
    res_ids:
        List of residue ids.
    atom_coord:
        Array of atom coordinates. Each residue is contained in a subarray.
    atom_symbols:
        List of atom symbols. Each residue is contained in a sublist.
    conn: np.array (N,cdim).
        connectivity matrix of the system.
    '''

    with open(pdbfile) as f:

        FoundAt = False
        FoundBond = False
        atom_idxs = []
        atom_names = []
        chain_idxs = []
        res_names = []
        res_ids = []
        atom_coords = []
        atom_symbols = []
        Ib1 = []
        Ib2 = []

        for line in f:

            # Read Atoms
            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':

                atom_idx = int(line[6:11])
                atom_name = line[11:16].strip()
                chain_idx = line[16].strip()
                res_name = line[17:20].strip()
                chain_idx = line[21].strip()
                res_id = int(line[22:26])
                coor = list(map(float, line[30:54].split()))
                atom_symbol = line[76:78].strip()

                # Special case for the first residue
                try:
                    res_ids[-1]
                except IndexError:
                    res_ids.append(res_id)
    
                # Special case for the first residue
                try:
                    res_names[-1]
                except IndexError:
                    res_names.append(res_name)
                
                # Check id: new residue or old one
                # if in new residue
                if res_id != res_ids[-1]:
                    res_ids.append(res_id)
                    res_names.append(res_name)

                    # save the old ones
                    atom_idxs.append(idxs)
                    atom_names.append(names)
                    atom_coords.append(np.array(coords))
                    atom_symbols.append(symbols)
                    
                    # initialize data for the new one
                    idxs = [atom_idx]
                    names = [atom_name]
                    coords = [coor]
                    symbols = [atom_symbol]

                # if still in the old residue
                else:
                
                    # save data for the new atom
                    try:
                        idxs.append(atom_idx)
                        names.append(atom_name)
                        coords.append(coor)
                        symbols.append(atom_symbol)
                
                    # unless no atom has been saved before
                    except:
                        idxs = [atom_idx]
                        names = [atom_name]
                        coords = [coor]
                        symbols = [atom_symbol]

            # Connectivity
            elif line[0:6] == 'CONECT':

                data = line.split()[1:]

                if len(data) > 1:
                    Ib1 += [int(data[0])]
                    tmpIb2 = np.zeros(4, dtype=int)
                    tmp1Ib2 = list(map(int, data[1:]))

                    for k in range(len(tmp1Ib2)):
                        tmpIb2[k] += tmp1Ib2[k]

                    Ib2 += [tmpIb2.tolist()]

        # save data for the last residue 
        atom_idxs.append(idxs)
        atom_names.append(names)
        atom_coords.append(np.array(coords))
        atom_symbols.append(symbols)

    atom_coords = np.array(atom_coords)

    # Build Connectivity Matrix
    NAtoms = sum([ len(i) for i in atom_idxs ])
    Conn = np.zeros((NAtoms, cdim), dtype=int)

    for i in range(len(Ib1)):
        Conn[Ib1[i] - 1] = Ib2[i]

    return atom_idxs, atom_names, res_names, res_ids, atom_coords, atom_symbols, Conn


def build_neigh_matrix(conn):
    '''
    Function to build a 1-2, 1-3, 1-4 neighbours matrix given the connectivity.

    Parameters
    ----------
    conn: np.array (N,M).
        connectivity matrix of the system.

    Returns
    -------
    neigh: np.array (N,N).
        Neighbours matrix.
    '''

    neigh = np.zeros((conn.shape[0],conn.shape[0]))

    for i in range(conn.shape[0]):
        neigh[i,i] = 1
        for j in range(conn.shape[1]):
            if conn[i,j] == 0:
                continue
            else:
                neigh[conn[i,j] - 1, i] = 2
                neigh[i, conn[i,j] - 1] = 2

    for i in range(neigh.shape[0]):
        for j in range(neigh.shape[0]):
            if neigh[i,j] == 2:
                for k in range(neigh.shape[0]):
                    if neigh[j,k] == 2 and neigh[i,k] != 1 \
                    and neigh[i,k] != 2:
                        neigh[i,k] = 3
                        neigh[k,i] = 3

    for i in range(neigh.shape[0]):
        for j in range(neigh.shape[0]):
            if neigh[i,j] == 3:
                for k in range(neigh.shape[0]):
                    if neigh[j,k] == 2 and neigh[i,k] != 1 \
                    and neigh[i,k] != 2 and neigh[i,k] != 3:
                        neigh[i,k] = 4
                        neigh[k,i] = 4

    return neigh


def build_R_matrix(coord, pol, a=1.7278):
    '''
    Function to fill in the relay matrix R (see JPCA 1998, 102, 2399).

    Parameters
    ----------
    coord: np.array (N,3).
        coordinates in au.
    pol: np.array (N).
        isotropic atomic polarisabilities in au.

    Returns
    -------
    R: np.array (3*N,3*N).
        Relay matrix.
    '''

    R = np.zeros((3 * coord.shape[0],3 * coord.shape[0]))

    for i in range(R.shape[0])[::3]:

        # Diagonal terms
        R[i:i+3,i:i+3] = np.eye(3) / pol[i // 3]

        # Non-diagonal terms, only upper half
        for j in range(i + 3, R.shape[0])[::3]:

            # Compute distance between sites i and j
            r = coord[j // 3] - coord[i // 3]
            rij = np.linalg.norm(r)

            # Compute screening distance between sites i and j
            s = a * ( (pol[i // 3] * pol[j // 3])**(1.0 / 6.0) )

            # Thole screening for linear charge distribution
            if rij <= s:
                v = rij / s
                s3 = 4 * v**3 - 3 * v**4
                s5 = v**4
            else:
                s3 = 1
                s5 = 1

            # Compute dipole field tensor between sites i and j
            Tij = ( s3 * np.eye(r.shape[0]) / rij**3 ) - ( s5 * 3 * np.outer(r, r) / rij**5 )
            R[i:i+3,j:j+3] = Tij

    # Symmetrise R
    R = R + R.T - np.diag(R.diagonal())

    return R


def reorganise_R_matrix(R):
    '''
    Function to reorganise a (3*N,3*N) matrix from an atomic-block diagonal
    (x_1, y_1, z_1...x_N, y_N, z_N) to a coordinate-block diagonal
    (x_1...x_N, y_1...y_N, z_1...z_N).

    Parameters
    ----------
    R: np.array (3*N,3*N).
        Relay matrix.

    Returns
    -------
    R: np.array (3*N,3*N).
        Relay matrix.
    '''

    for i in range(3):
        try:
            idxs = np.concatenate((idxs, np.arange(i, R.shape[0])[::3]))
        except:
            idxs = np.arange(i, R.shape[0])[::3]

    R = R[idxs][:,idxs]

    return R


def kabsch(struct1, struct2):
    '''Returns the RMSD calculated with Kabsch's algorithm.'''

    # Modify structures to get rid of the atomic symbol or number and convert
    # to np.array
    struct1 = np.array([ [atom[0], atom[1], atom[2]] for atom in struct1 ])    
    struct2 = np.array([ [atom[0], atom[1], atom[2]] for atom in struct2 ])    

    # check for consistency in number of atoms
    assert len(struct1) == len(struct2)
    L = len(struct1)
    assert L > 0

    # Center the two fragments to their center of coordinates
    com1 = np.sum(struct1, axis=0) / float(L)
    com2 = np.sum(struct2, axis=0) / float(L)
    struct1 -= com1
    struct2 -= com2

    # Initial residual, see Kabsch.
    E0 = np.sum(np.sum(struct1 * struct1, axis=0), axis=0) + \
         np.sum(np.sum(struct2 * struct2, axis=0), axis=0)

    # This beautiful step provides the answer. V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    V, S, Wt = np.linalg.svd(np.dot(np.transpose(struct2), struct1))

    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation. V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    # The rotation matrix U is simply V*Wt
    U = np.dot(V, Wt)
 
    # rotate and translate the molecule
    struct2 = np.dot((struct2), U)
    struct2 = struct2 + com1

    return struct2, U


def banner(text=None, ch='=', length=78):
    '''
    Return a banner line centering the given text.
    
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
    '''
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
