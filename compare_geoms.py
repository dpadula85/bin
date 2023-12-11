#!/usr/bin/env python

import os
import csv
import sys
import numpy as np
import pandas as pd
import networkx as nx
import argparse as arg
import MDAnalysis as mda
from MDAnalysis.lib.util import unique_rows
from rdkit.Chem.rdmolfiles import MolFromXYZFile
from rdkit.Chem.AllChem import GetBestRMS, AlignMol

import warnings
warnings.filterwarnings("ignore")

au2ang = 0.5291771


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument(
            '-r',
            '--ref',
            type=str,
            required=True,
            dest='RefGeom',
            help='''Reference Geometry.'''
        )

    inp.add_argument(
            '-g',
            '--geom',
            type=str,
            required=True,
            dest='Geom',
            help='''Geometry to be compared to the reference.'''
        )

    #
    # Outut Options
    #
    out = parser.add_argument_group("Outut Data")

    out.add_argument(
            '-o',
            '--out',
            type=str,
            default=None,
            dest='OutFile',
            help='''Output file prefix.'''
        )

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


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


def rmse(x, y):
    '''
    Function to compute the Root Mean Square Error between two sets of points.

    Parameters
    ----------
    x: np.array (N,M).
        set of points.
    y: np.array (N,M).
        set of points.

    Returns
    -------
    err: float.
        Root Mean Squared Error.
    '''

    err = np.sqrt(np.mean( (x - y)**2))

    return err


def get_sp2(u):
    '''
    Function to get indices of sp2 and sp3 atoms.

    Parameters
    ----------
    u: object.
        MDAnalysis Universe to be cropped.

    Returns
    -------
    sp2: np.ndarray.
        Indices of sp2 atoms.
    sp3: np.ndarray.
        Indices of sp3 atoms.
    '''

    # Get bonds
    try:
        bds = u.bonds.to_indices()
    except:
        u.guess_bonds()
        bds = u.bonds.to_indices()

    # Get connectivity, use -1 as a placeholder for empty valence
    conn = np.ones((len(u.atoms), 4)) * -1
    for bond in bds:
        at1, at2 = bond
        for j in np.arange(conn[at1].shape[0]):
            if conn[at1,j] == -1:
                conn[at1,j] = at2
                break

        for j in np.arange(conn[at2].shape[0]):
            if conn[at2,j] == -1:
                conn[at2,j] = at1
                break

    # Get heavy atoms
    heavy = np.where(u.atoms.types != "H")[0]

    # Get sp3 atoms
    sp3 = np.where(np.all(conn > -1, axis=1))[0]
    allcheck = sp3.copy()

    # Check all sp3 atoms
    keep = []
    delete = []
    for satat in allcheck:

        # check connectivity
        iconn = conn[satat]

        # filter H out from connected
        iconnheavy = iconn[np.in1d(iconn, heavy)]

        delete.append(satat)
        delete.extend(iconn[~np.in1d(iconn, heavy)])

    # Convert to int arrays
    keep = np.asarray(keep).astype(int)
    delete = np.asarray(delete).astype(int)

    # Get non sp3 atoms
    unsat = ~np.all(conn > -1, axis=1)

    # Set which saturated atoms to keep or delete
    unsat[keep] = True
    unsat[delete] = False
    tokeep = np.where(unsat)[0]
    sp2 = np.intersect1d(tokeep, heavy)

    return sp2, sp3


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


def bond(A, B):
    '''
    Function to compute the distance defined by points A and B.

    Parameters
    ----------
    A: np.array (N).
        First point.
    B: np.array (N).
        Second point.

    Returns
    -------
    ab: float.
        Distance.
    '''

    ab = np.linalg.norm(B - A)

    return ab


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
    dihedral = np.degrees(np.arctan2(y, x))

    return dihedral


def findPaths(G, u, n):
    '''
    Function to get n-order paths in graph G, starting from node u.

    Parameters
    ----------
    G: NetworkX Graph object.
        Graph for path search.
    u: int.
        Index of the starting node for the path search.
    n: int.
        Order of the paths to be found, that is number of bonds separating
        atoms involved in the interaction (1 for bonds, 2 for angles, 3 for
        dihedrals or 1,4 interactions, 4 for 1,5 interactions etc).

    Returns
    -------
    paths: list.
        List of n-order paths found in Graph G starting from node u.
    '''

    if n == 0:
        return [[u]]
    paths = [
        [u] + path for neighbor in G.neighbors(u) for path 
        in findPaths(G,neighbor,n-1) if u not in path
    ]

    return paths


def flip_bigger(arr, col1=0, col2=-1):
    '''
    Function to flip columns of array arr based on values in columns
    identified by index col1 and col2.

    Parameters
    ----------
    arr: np.ndarray (shape: (N, M)).
        Array whose columns are to be flipped.
    col1: int.
        Index of the first column.
    col2: int.
        Index of the second column.

    Returns
    -------
    arr: np.ndarray (shape: (N, M)).
        Array with flipped columns (no duplicate rows).
    '''

    arr = arr.copy()
    flipped = arr[:,col1] > arr[:,col2]
    arr[flipped] = np.fliplr(arr[flipped])
    arr = unique_rows(arr)

    return arr


def myround(x, base=5):
    '''
    Function to round a number x to the closest multiple of an arbitrary base.

    Parameters
    ----------
    x: float.
        number to be rounded.
    base: int or float.
        Base for the closest multiple.

    Returns
    -------
    rounded: int or float.
        Rounded number to the closest multiple of base.
    '''

    rounded = base * round(x / base)

    return rounded


def select_dihedrals(u, diheds, sp2, rings):
    '''
    Function to separate stiff and flexible dihedrals from the complete set
    based on the atoms making up a specific coordinate belonging to the same
    ring. Additionally, it selects the stiff dihedrals to keep (the ones at
    0 degrees), discarding the redundant set (the ones at 180 degrees).

    Parameters
    ----------
    u: object.
        MDAnalysis Universe, to get coordinates.
    diheds: np.ndarray (shape: (N, 4))
        Array of bonds, in terms of indices of involved atoms.
    sp2: np.ndarray.
        Indices of sp2 atoms.
    eq_hs: dict.
        Dictionary of equivalent H atoms.
    rings: dict.
        Dictionary of rings, divided in core atoms, and core plus substituents
        (full).

    Returns
    -------
    stiff: np.ndarray (shape: (N, 4)).
        Array of stiff proper dihedrals, in terms of indices of involved atoms.
    flex: np.ndarray (shape: (N, 4)).
        Array of flexible proper dihedrals, in terms of indices of involved
        atoms.
    '''

    # Divide dihedrals into stiff and flexible
    # if the two central atoms belong to the same ring, then stiff
    # otherwise flexible. Only need to check proper dihedrals
    flex = []
    stiff = []
    for dihed in diheds:
        ati, atj, atk, atl = dihed
        for ring in rings['core']:
            if atj in ring and atk in ring and atj in sp2 and atk in sp2:
                stiff.append(dihed)
                break

        try:
            check = (dihed == np.asarray(stiff)).all(axis=1).any()
        except:
            check = False

        if not check:
            flex.append(dihed)

    flex = np.asarray(flex)
    stiff = np.asarray(stiff)

    return stiff, flex


def list_intcoords(coordfile):
    '''
    Function to obtain a set of internal coordinates from a geometry file.

    Parameters
    ----------
    coordfile: str.
        Name of the geometry file for which internal coordinates are required.

    Returns
    -------
    bds: np.ndarray.
        Array of bonds, in terms of indices of involved atoms.
    angles: np.ndarray.
        Array of angles, in terms of indices of involved atoms.
    stiff: np.ndarray (shape: (N, 4)).
        Array of stiff proper dihedrals, in terms of indices of involved atoms.
    impdiheds: np.ndarray.
        Array of improper dihedrals, in terms of indices of involved atoms.
    flex: np.ndarray (shape: (N, 4)).
        Array of flexible proper dihedrals, in terms of indices of involved
        atoms.
    rings: dict.
        Dictionary of rings, divided in core atoms, and core plus substituents
        (full).
    '''

    u = mda.Universe(coordfile, guess_bonds=True)

    # Get bonds
    bds = u.bonds.to_indices()

    # Build a graph representation to easily find all other
    # internal coordinates
    nodes = np.arange(bds.min(), bds.max() + 1)
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(bds)

    # Get rings
    rings = {
        'core' : sorted(nx.cycle_basis(G))
    }

    # Get angles
    angles = []
    for node in nodes:
        angs = findPaths(G, node, 2)
        angles.extend(angs)

    angles = np.asarray(angles)

    # make first entry always less than third
    angles = flip_bigger(angles)

    # order by central atom and then by first
    angles = angles[angles[:,0].argsort()]
    angles = angles[angles[:,1].argsort(kind='mergesort')]

    # Get proper dihedrals
    diheds = []
    for node in nodes:
        dihs = findPaths(G, node, 3)
        diheds.extend(dihs)

    diheds = np.asarray(diheds)

    # make first entry always less than last
    diheds = flip_bigger(diheds, col1=1, col2=2)

    # order by central pair and then by first
    diheds = diheds[diheds[:,0].argsort()]
    diheds = diheds[diheds[:,1].argsort(kind='mergesort')]
    diheds = diheds[diheds[:,2].argsort(kind='mergesort')]

    # Get sp2 and sp3 heavy atoms to select stiff and look for impropers
    sp2, sp3 = get_sp2(u)

    # Separate stiff from flexible and select which stiff to keep
    stiff, flex = select_dihedrals(u, diheds, sp2, rings)

    # Get improper dihedrals
    impdiheds = []
    for idx in sp2:
        try:
            neighs = [ i for i in nx.all_neighbors(G, idx) ]
        except:
            neighs = []
        if len(neighs) > 2:
            impdiheds.append([ idx ] + neighs)

    impdiheds = np.asarray(impdiheds)

    # order by first atom
    # exception to deal with no impropers
    try:
        impdiheds = impdiheds[impdiheds[:,0].argsort()]
    except:
        pass

    return bds, angles, stiff, impdiheds, flex, rings


def get_int_coords(coordfile, **kwargs):
    '''
    Function to obtain the values of internal coordinates.

    Parameters
    ----------
    coordfile: str.
        Name of the geometry file for which internal coordinates are required.
    kwargs: dict.
        Dictionary of internal coordinates to evaluate.

    Returns
    -------
    bdsv: np.ndarray.
        Array of bond lengths.
    angsv: np.ndarray.
        Array of bond angles.
    stiffv: np.ndarray.
        Array of stiff proper dihedrals.
    impdihedsv: np.ndarray.
        Array of improper dihedrals.
    flexv: np.ndarray.
        Array of flexible proper dihedrals.
    '''

    u = mda.Universe(coordfile)
    bds, angles, stiff, impdiheds, flex, rings = list_intcoords(coordfile)
    bds = kwargs.pop("bds", bds)
    angles = kwargs.pop("angles", angles)
    stiff = kwargs.pop("stiff", stiff)
    impdiheds = kwargs.pop("impdiheds", impdiheds)
    flex = kwargs.pop("flex", flex)

    bds_x = u.atoms.positions[bds]
    bdsv = np.linalg.norm(bds_x[:,1,:] - bds_x[:,0,:], axis=1)

    angs_x = u.atoms.positions[angles]
    angsv = np.asarray([ angle(*i) for i in angs_x ])

    try:
        stiff_x = u.atoms.positions[stiff]
        stiffv = np.asarray([ dihedral(*i) for i in stiff_x ])
    except IndexError:
        stiffv = None

    try:
        impdiheds_x = u.atoms.positions[impdiheds]
        impdihedsv = np.asarray([ dihedral(*i) for i in impdiheds_x ])
    except IndexError:
        impdihedsv = None

    try:
        flex_x = u.atoms.positions[flex]
        flexv = np.asarray([ dihedral(*i) for i in flex_x ])
    except IndexError:
        flexv = None

    return bdsv, angsv, stiffv, impdihedsv, flexv


def compare_geoms(coordfile1, coordfile2):
    '''
    Function to compare two geometries in terms of internal coordinates.

    Parameters
    ----------
    coordfile1: str.
        Name of the first geometry file.
    coordfile2: str.
        Name of the second geometry file.

    Returns
    -------
    bds_rmse: float.
        RMSE between bond lengths (in Angstroem).
    angs_rmse: float.
        RMSE between bond angles (in degrees).
    stiff_rmse: float.
        RMSE between stiff dihedrals (in degrees).
    impdiheds_rmse: float.
        RMSE between improper dihedrals (in degrees).
    flex_rmse: float.
        RMSE between flexible dihedrals (in degrees).
    tot_rmse: float.
        RMSE between atomic coordinates (in Angstroem).
    '''

    bds, angles, stiff, impdiheds, flex, rings = list_intcoords(coordfile1)
    int_coords = {
        "bds" : bds,
        "angles" : angles,
        "stiff" : stiff,
        "impdiheds" : impdiheds,
        "flex" : flex,
    }
    bdsv1, angsv1, stiffv1, impdihedsv1, flexv1 = get_int_coords(coordfile1)
    bdsv2, angsv2, stiffv2, impdihedsv2, flexv2 = get_int_coords(coordfile2, **int_coords)

    bds_rmse = rmse(bdsv1, bdsv2)
    angs_rmse = rmse(angsv1, angsv2)
    try:
        stiff_rmse = rmse(stiffv1, stiffv2)
    except TypeError:
        stiff_rmse = None

    try:
        impdiheds_rmse = rmse(impdihedsv1, impdihedsv2)
    except TypeError:
        impdiheds_rmse = None

    try:
        flex_rmse = rmse(flexv1, flexv2)
    except TypeError:
        flex_rmse = None

    mol1 = MolFromXYZFile(coordfile1)
    mol2 = MolFromXYZFile(coordfile2)
    tot_rmse = AlignMol(mol1, mol2)

    return bds_rmse, angs_rmse, stiff_rmse, impdiheds_rmse, flex_rmse, tot_rmse


if __name__ == '__main__':

    opts = options()

    data = compare_geoms(opts['RefGeom'], opts['Geom'])
    data = np.asarray(list(data)).reshape(1, -1)

    cols = [
        "bonds",
        "angles",
        "stiff",
        "impropers",
        "flexible",
        "total"
    ]

    df = pd.DataFrame(data=data, columns=cols)

    if opts['OutFile'] is not None:
        df.to_csv(f"{opts['OutFile']}.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)
    else:
        print(df)
