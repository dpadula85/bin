#!/usr/bin/env python

import csv
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import argparse as arg
import MDAnalysis as mda
from scipy.spatial.distance import cdist


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument('--top', type=str, dest='TopFile', required=True,
                     help='''Topology File.''')

    inp.add_argument('-ds', '--dsel', type=str, dest='DSel', required=True,
                     help='''Donor Selection.''')

    inp.add_argument('-as', '--asel', type=str, dest='ASel', required=True,
                     help='''Acceptor Selection.''')

    inp.add_argument('--trj', type=str, dest='TrjFile', default=None,
                     help='''Trajectory File.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-p', '--pre', default="dimer_analysis", type=str,
                     dest='OutPre', help='''Output File Prefix.''')

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


def euler_angles_from_matrix(R):
    '''
    See http://www.close-range.com/docs/Computing_Euler_angles_from_a_rotation_matrix.pdf

    Parameters
    ----------
    R: np.array (3,3).
        Rotation matrix.

    Returns
    -------
    psi: float.
        Roll angle (in degrees).
    theta: float.
        Pitch angle (in degrees).
    phi: float.
        Yaw angle (in degrees).
    '''

    phi = 0.0
    if np.isclose(R[2,0], -1.0):
        theta = np.pi / 2.0
        psi = np.arctan2(R[0,1], R[0,2])

    elif np.isclose(R[2,0], 1.0):
        theta = -np.pi / 2.0
        psi = np.arctan2(-R[0,1], -R[0,2])

    else:
        theta = -np.arcsin(R[2,0])
        cos_theta = np.cos(theta)
        psi = np.arctan2(R[2,1] / cos_theta, R[2,2] / cos_theta)
        phi = np.arctan2(R[1,0]/cos_theta, R[0,0] / cos_theta)

    return np.degrees(psi), np.degrees(theta), np.degrees(phi)


def find_all_rings(sel):
    '''
    Function to find rings in a structure.

    Parameters
    ----------
    sel: MDAnalysis AtomGroup.
        Group selection to look for rings.

    Returns
    -------
    rings: list of lists.
        list of indices of rings in the input selection.
    '''

    # Handle atoms as Graph
    g = nx.Graph()
    g.add_edges_from(sel.bonds.to_indices())

    # cycle_basis gives all rings
    cycles = nx.cycle_basis(g)

    # select rings only
    ring_members = set(itertools.chain(*cycles))
    idxs = "index %s" % ' '.join(map(str, ring_members))
    ring_ag = sel.select_atoms(idxs)

    # make a new graph, this time of all rings
    g = nx.Graph()
    for cycle in cycles:
        nx.function.add_cycle(g, cycle)

    # Find connected parts of the graph (manages fused rings)
    rings = [ list(i) for i in nx.connected_components(g) ]

    return rings


def analyse_dimer(donor, accpt):
    '''
    Function to compute geometrical quantities for a D/A pair, based on the
    principal axes of inertia.

    Parameters
    ----------
    donor: MDAnalysis AtomGroup.
        Donor atoms
    accpt: MDAnalysis AtomGroup.
        Acceptor atoms.

    Returns
    -------
    rcoms: float.
        Distance between the centres of mass of the two groups (in Angstroem).
    rpi: float.
        Pi-stacking distance between the pi stacked rings (in Angstroem).
    psi: float.
        Roll angle (in degrees).
    theta: float.
        Pitch angle (in degrees).
    phi: float.
        Yaw angle (in degrees).
    '''

    donor.guess_bonds()
    accpt.guess_bonds()

    # donor COM and principal axes
    # Invert order of principal axes so to have the long molecular axis first,
    # short second and the orthogonal one last, and check righthandedness.
    dcom = donor.center_of_mass()
    dpa = donor.principal_axes()[::-1].T
    if np.linalg.det(dpa) < 0:
        dpa[:,-1] = -dpa[:,-1]

    # accpt COM and principal axes
    # Invert order of principal axes so to have the long molecular axis first,
    # short second and the orthogonal one last, and check righthandedness.
    acom = accpt.center_of_mass()
    apa = accpt.principal_axes()[::-1].T
    if np.linalg.det(apa) < 0:
        apa[:,-1] = -apa[:,-1]

    # Define rotation matrix between the two sets of principal axes.
    # M transforms from the dpa frame to the apa frame.
    M = np.dot(apa, np.linalg.inv(dpa))

    # Get roll, pitch, yaw angles associated to the transformation matrix
    # Pitch will be recomputed once pi-stacking sites will be identified
    psi, theta, phi = euler_angles_from_matrix(M)

    # Transform between -90 and 90
    psi = np.degrees(( np.radians(psi) + np.pi / 2 ) % np.pi - np.pi / 2)
    phi = np.degrees(( np.radians(phi) + np.pi / 2 ) % np.pi - np.pi / 2)

    # Distance between coms
    r = acom - dcom
    rcoms = np.linalg.norm(r)

    # Find all rings in donor and get their COMs
    drings = find_all_rings(donor)
    dsels = [ "index %s" % ' '.join(map(str, x)) for x in drings ]
    drings_coms = [ donor.select_atoms(sel).center_of_mass() for sel in dsels ]
    drings_coms = np.array(drings_coms)

    # Find all rings in accpt and get their COMs
    arings = find_all_rings(accpt)
    asels = [ "index %s" % ' '.join(map(str, x)) for x in arings ]
    arings_coms = [ accpt.select_atoms(sel).center_of_mass() for sel in asels ]
    arings_coms = np.array(arings_coms)

    # Find two closest rings
    D = cdist(drings_coms, arings_coms)
    didx, aidx = np.unravel_index(D.argmin(), D.shape)

    # Get their coordinates
    dring = donor.select_atoms(dsels[didx])
    aring = accpt.select_atoms(asels[aidx])

    # Pi-stacking distance as minimum distance between atoms in the two rings
    rpi = cdist(dring.atoms.positions, aring.atoms.positions).min()

    # Fit a plane to each ring and get the unit vectors
    ndring = lstsq_fit(dring.atoms.positions)
    ndring /= np.linalg.norm(ndring)
    naring = lstsq_fit(aring.atoms.positions)
    naring /= np.linalg.norm(naring)

    # Compute angle between unit vectors normal to the planes
    theta = v1v2_angle(ndring, naring)
    theta = np.degrees(( np.radians(theta) + np.pi / 2 ) % np.pi - np.pi / 2)

    return rcoms, rpi, psi, theta, phi


def main(**Opts):

    if Opts["TrjFile"]:
        u = mda.Universe(Opts["TopFile"], Opts["TrjFile"])
    else:
        u = mda.Universe(Opts["TopFile"])

    donor = u.select_atoms(Opts["DSel"])
    accpt = u.select_atoms(Opts["ASel"])
    data = []
    for ts in u.trajectory:
        t = ts.time
        r, rpi, roll, pitch, yaw = analyse_dimer(donor, accpt)
        snapdata = np.array([ t, r, rpi, roll, pitch, yaw ])
        data.append(snapdata)

    data = np.array(data)
    df = pd.DataFrame({
        "Time / ns": data[:,0] / 1000.0,
        "r / A" : data[:,1],
        "rpi / A" : data[:,2],
        "roll / deg" : data[:,3],
        "pitch / deg" : data[:,4],
        "yaw / deg" : data[:,5],
        })

    df.to_csv("%s.csv" % Opts["OutPre"], quoting=csv.QUOTE_NONNUMERIC, index=False)

    return df


if __name__ == '__main__':
    Opts = options()
    main(**Opts)
