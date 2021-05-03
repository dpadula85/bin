#!/usr/bin/env python

import csv
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import argparse as arg
import MDAnalysis as mda


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
    rings: list of MDAnalysis AtomGroups.
        list of rings in the input selection.
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
    rings_idxs = [ list(i) for i in nx.connected_components(g) ]

    rings = [ sel.select_atoms("index %s" % ' '.join(map(str, i)))
              for i in rings_idxs if len(i) > 4 ]

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
    rnorm: float.
        Distance between the centres of mass of the two groups (in Angstroem).
    rd: float.
        Longitudinal displacement between the two groups (in Angstroem).
    rs: float.
        Side displacement between the two groups (in Angstroem).
    rpi: float.
        Pi-Stacking distance between the two groups (in Angstroem).
    psi: float.
        Roll angle (in degrees).
    theta: float.
        Pitch angle (in degrees).
    phi: float.
        Yaw angle (in degrees).
    '''

    # Donor quantities
    donor.guess_bonds()
    drings = find_all_rings(donor)
    didxs = flatten([ x.indices.tolist() for x in drings ])
    didxs = "index %s" % ' '.join(map(str, didxs))
    donor = donor.select_atoms(didxs)

    # COM
    dcom = donor.center_of_mass()

    # Invert order of principal axes so to have the long molecular axis first,
    # short second and the orthogonal one last, and check righthandedness.
    dpa = donor.principal_axes()[::-1].T
    if np.linalg.det(dpa) < 0:
        dpa[:,-1] = -dpa[:,-1]

    # Acceptor quantities
    accpt.guess_bonds()
    arings = find_all_rings(accpt)
    aidxs = flatten([ x.indices.tolist() for x in arings ])
    aidxs = "index %s" % ' '.join(map(str, aidxs))
    accpt = accpt.select_atoms(aidxs)

    # COM
    acom = accpt.center_of_mass()

    # Invert order of principal axes so to have the long molecular axis first,
    # short second and the orthogonal one last, and check righthandedness.
    apa = accpt.principal_axes()[::-1].T
    if np.linalg.det(apa) < 0:
        apa[:,-1] = -apa[:,-1]

    # Define rotation matrix between the two sets of principal axes.
    # M transforms from the dpa frame to the apa frame.
    M = np.dot(apa, np.linalg.inv(dpa))

    # Get roll, pitch, yaw angles associated to the transformation matrix
    psi, theta, phi = euler_angles_from_matrix(M)

    # Transform between -90 and 90
    psi = np.degrees(( np.radians(psi) + np.pi / 2 ) % np.pi - np.pi / 2)
    theta = np.degrees(( np.radians(theta) + np.pi / 2 ) % np.pi - np.pi / 2)
    phi = np.degrees(( np.radians(phi) + np.pi / 2 ) % np.pi - np.pi / 2)

    # Distance between coms
    r = acom - dcom
    rnorm = np.linalg.norm(r)

    # Project distance components onto the accpt principal axes
    # by convention, because acceptors are smaller, thus they
    # are better described by a plane
    rd, rs, rpi = np.abs(np.dot(apa, r))

    return rnorm, rd, rs, rpi, psi, theta, phi


def main(**Opts):

    if Opts["TrjFile"]:
        u = mda.Universe(Opts["TopFile"], Opts["TrjFile"])
    else:
        u = mda.Universe(Opts["TopFile"])

    data = []
    for ts in u.trajectory:
        t = ts.time
        donor = u.select_atoms(Opts["DSel"])
        # donor.guess_bonds()
        accpt = u.select_atoms(Opts["ASel"])
        # accpt.guess_bonds()
        r, rd, rs, rpi, roll, pitch, yaw = analyse_dimer(donor, accpt)
        snapdata = np.array([ t, r, rd, rs, rpi, roll, pitch, yaw ])
        data.append(snapdata)

    data = np.array(data)
    df = pd.DataFrame({
        "Time / ns": data[:,0] / 1000.0,
        "r / A" : data[:,1],
        "rd / A" : data[:,2],
        "rs / A" : data[:,3],
        "rpi / A" : data[:,4],
        "roll / deg" : data[:,5],
        "pitch / deg" : data[:,6],
        "yaw / deg" : data[:,7],
        })

    df.to_csv("%s.csv" % Opts["OutPre"], quoting=csv.QUOTE_NONNUMERIC, index=False)

    return df


if __name__ == '__main__':
    Opts = options()
    main(**Opts)
