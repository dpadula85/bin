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
    Function to obtain the connectivity matrix of a molecule from the
    coordinates and the atomic radii.

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

    # Add fictitious distance on the diagonal not to find the atom itself
    D = D + np.diag(D.diagonal() + np.inf)

    # Compute radii sum matrix
    R = radii.reshape(1, -1) / 2 + radii.reshape(-1, 1) / 2

    # Compare them to figure out connectivity
    idxs = D <= R

    # For each line get True idxs
    connidxs = [ np.where(i == True)[0].tolist() for i in idxs ]
    pad = len(max(connidxs, key=len))
    M = np.array([ i + [ -1 ] * (pad - len(i)) for i in connidxs ])

    return M


def rot(axis, theta):
    '''
    Returns the rotation matrix for the anticlockwise rotation about an
    arbitrary axis by theta according to Rodrigues' formula.

    Parameters
    ----------
    axis: np.array (N).
        Unit vector describing the rotation axis.
    theta: float.
        Angle of rotation (in degrees).

    Returns
    -------
    R: np.array (N,N).
        Rotation matrix.
    '''

    axis = axis / np.linalg.norm(axis)
    theta = np.radians(theta)
    I = np.eye(axis.shape[0])

    # Define axis' cross-product matrix
    K = np.cross(I, axis)

    R = I + np.sin(theta) * K + (1 - np.cos(theta)) * np.linalg.matrix_power(K, 2)

    return R


def v1v2_angle(v1, v2):
    '''
    Function to compute the angle between two vectors.

    Parameters
    ----------
    v1: np.array (N).
        First vector.
    v2: np.array (N).
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


def rot_vecs(v1, v2):
    '''
    Returns the rotation matrix that transforms v1 into v2.

    Parameters
    ----------
    v1: np.array (N).
        First vector.
    v2: np.array (N).
        Second vector.

    Returns
    -------
    R: np.array (N,N).
        Rotation matrix.
    '''

    # Find rotation vector orthogonal to the plane spanned by v1 and v2
    axis = np.cross(v1, v2)
    axis = axis / np.linalg.norm(axis)

    # Find the angle of rotation
    theta = v1v2_angle(v1, v2)

    # Obtain the rotation matrix with Rodrigues' formula
    R = rot(axis, theta)

    return R


def angle(A, B, C):
    '''
    Function to compute the angle defined by points A, B, and C.

    Parameters
    ----------
    A: np.array (N).
        First point.
    B: np.array (N).
        Second point, vertex of the angle.
    C: np.array (N).
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
    A: np.array (N).
        First point.
    B: np.array (N).
        Second point.
    C: np.array (N).
        Third point.
    D: np.array (N).
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


def lstsq_fit(pts):
    '''
    Function to fit a set of points with least squares. The geometrical objects
    involved depend on the number of dimensions (columns) of pts.

    Parameters
    ----------
    pts: np.array (N,M).
        coordinates.

    Returns
    -------
    coeffs: np.array (M).
        coefficients of the least squares fit.
    '''

    A = np.c_[ pts[:,:-1], np.ones(pts.shape[0]) ]
    B = pts[:,-1]

    coeffs, res, rank, singular_values = np.linalg.lstsq(A, B, rcond=None)

    return coeffs


def make_grid(**kwargs):
    '''
    Function to create a grid of points centred at the origin according to
    the basis vectors. Optional kwargs can control both the reference frame
    and whether the grid will be generated in space, onto a plane, or along a
    line.

    Parameters
    ----------
    ref: np.array (3,3).
        basis vectors.
    origin: np.array (3).
        coordinated of the origin.
    xu, yu, zu: float.
        maximum coefficient for each basis vector.
    xl, yl, zl: float.
        minimum coefficient for each basis vector.
    nx, ny, nz: int.
        number of points along each direction.

    Returns
    -------
    grid: np.array (M,3).
        grid of points.
    '''

    # Assign default reference system and origin to the cartesian ref frame
    ref = kwargs.pop("ref", np.eye(3))
    origin = kwargs.pop("origin", np.zeros(3))

    # Define options for the grid
    xu = kwargs.pop("xu", 5)
    yu = kwargs.pop("yu", 5)
    zu = kwargs.pop("zu", 5)
    xl = kwargs.pop("xl", -xu)
    yl = kwargs.pop("yl", -yu)
    zl = kwargs.pop("zl", -zu)
    nx = kwargs.pop("nx", np.abs(xl) + np.abs(xu) + 1)
    ny = kwargs.pop("ny", np.abs(yl) + np.abs(yu) + 1)
    nz = kwargs.pop("nz", np.abs(zl) + np.abs(zu) + 1)

    # Define spacings along each basis vector
    i = np.linspace(xl, xu, nx)
    j = np.linspace(yl, yu, ny)
    k = np.linspace(zl, zu, nz)

    # We should do
    # for p in i:
    #     for q in j:
    #         for r in k:
    #             gridpoint = origin + p * e1 + q * e2 + r * e3
    #
    # where e1, e2, e3 are basis vectors stored as columns of ref

    # Here is a vectorised version of the nested for loop
    # Make grid of displacements along each basis vector
    g = np.meshgrid(i, j, k)

    # Convert to a more natural format, one column for each basis vector
    grid = np.vstack(list(map(np.ravel, g))).T

    # Transform each grid point to local basis and translate
    grid = np.dot(ref, grid.T).T + origin

    return grid


def proj(pts, coeffs):
    '''
    Function to compute projections of a set of points onto a plane.

    Parameters
    ----------
    pts: np.array (N,M).
        set of points
    coeffs: np.array (M).
        normal vector describing the plane.

    Returns
    -------
    prjs: np.array (N,M).
        set of points projected onto the plane.
    '''

    # For each point p in pts we should do
    # prj = p - np.dot(p, u) * u
    # where u is the normal unit vector

    # Compute normal vector
    u = coeffs / np.linalg.norm(coeffs)

    # Here is a vectorised version of the projection formula
    # Compute elementwise dot products between pts and u.
    # Deal with the case of a single point as an exception.
    try:
        dotprods = np.einsum('ij,ij->i', pts, u.reshape(1,-1))
    except:
        dotprods = np.einsum('ij,ij->i', pts.reshape(-1,u.shape[0]), u.reshape(1,-1))

    # Repeat dot products and unit vector for elementwise mult to be
    # subtracted from originary coordinates
    dotprodsmat = np.repeat(dotprods.reshape(-1,1), u.shape[0], axis=1)
    umat = np.repeat(u.reshape(1,-1), dotprods.shape[0], axis=0)

    # Subtract the components along the normal to the plane from the
    # originary coordinates
    prjs = pts - dotprodsmat * umat

    return prjs


def make_plane_basis(plane):
    '''
    Function to define a local reference frame on a plane defined by its
    normal vector.

    Parameters
    ----------
    plane: np.array (N).
        normal vector describing the plane.

    Returns
    -------
    ref: np.array (N,N).
        local basis.
    '''

    # Get missing coeff
    u = plane / np.linalg.norm(plane)
    p = np.zeros(3)
    p[-1] = plane[-1]
    d = np.dot(-p, u)

    # Define two points on the plane
    p1 = np.array([0, 0, d / u[2]])
    p2 = np.array([0, d / u[1], 0])
    e1 = p2 - p1
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(u, e1)

    # Define the local reference frame basis and check its righthandedness
    ref = np.c_[ e1, e2, u ]
    det = np.linalg.det(ref)
    if det < 0.0:
        ref[:,0] = -ref[:,0]

    return ref


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
    acf_t = np.array(list(map(r, t))) / c0

    return acf_t


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
            # These damping factors avoid polarisation catastrophe
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
