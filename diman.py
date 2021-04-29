#!/usr/bin/env python

import csv
import numpy as np
import pandas as pd
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
    rs: float.
        Pi-Stacking distance between the two groups (in Angstroem).
    ld: float.
        Longitudinal displacement between the two groups (in Angstroem).
    sd: float.
        Side displacement between the two groups (in Angstroem).
    alpha: float.
        Yaw angle (between long molecular axes, in degrees).
    beta: float.
        Roll angle (between short molecular axes, in degrees).
    gamma: float.
        Pitch angle (between vectors orthogonal to the molecular plane, in
        degrees).
    '''

    # Donor quantities
    # Invert order of principal axes so to have the long molecular axis first,
    # short second and the orthogonal one last
    dcom = donor.center_of_mass()
    dpa = donor.principal_axes()[::-1].T

    # Acceptor quantities
    # Invert order of principal axes so to have the long molecular axis first,
    # short second and the orthogonal one last
    acom = accpt.center_of_mass()
    apa = accpt.principal_axes()[::-1].T

    # Distance between coms
    r = acom - dcom
    rnorm = np.linalg.norm(r)

    # Yaw angle, angle between long axes (helicity)
    # reported between -90 and 90
    alpha = v1v2_angle(dpa[:,0], apa[:,0])
    alpha = np.degrees(( np.radians(alpha) + np.pi / 2 ) % np.pi - np.pi / 2)

    # Roll angle, angle between short axes
    # reported between -90 and 90
    beta = v1v2_angle(dpa[:,1], apa[:,1])
    beta = np.degrees(( np.radians(beta) + np.pi / 2 ) % np.pi - np.pi / 2)

    # Pitch angle, angle between axes orthogonal to molecular plane
    # reported between -90 and 90
    gamma = v1v2_angle(dpa[:,2], apa[:,2])
    gamma = np.degrees(( np.radians(gamma) + np.pi / 2 ) % np.pi - np.pi / 2)

    # Stacking distance, r projection on axis orthogonal to molecular plane
    rs1 = np.abs(np.dot(r, dpa[:,2]))
    rs2 = np.abs(np.dot(r, apa[:,2]))
    rs = np.min([ rs1, rs2 ])

    # Longitudinal displacement, r projection on long molecular axis.
    ld1 = np.abs(np.dot(r, dpa[:,0]))
    ld2 = np.abs(np.dot(r, apa[:,0]))
    ld = np.min([ ld1, ld2 ])

    # Side displacement, r projection on short molecular axis.
    sd1 = np.abs(np.dot(r, dpa[:,0]))
    sd2 = np.abs(np.dot(r, apa[:,0]))
    sd = np.min([ sd1, sd2 ])

    return rnorm, rs, ld, sd, alpha, beta, gamma


def main():

    Opts = options()
    if Opts["TrjFile"]:
        u = mda.Universe(Opts["TopFile"], Opts["TrjFile"])
    else:
        u = mda.Universe(Opts["TopFile"])

    data = []
    for ts in u.trajectory:
        t = ts.time
        donor = u.select_atoms(Opts["DSel"])
        accpt = u.select_atoms(Opts["ASel"])
        r, rs, ld, sd, alpha, beta, gamma = analyse_dimer(donor, accpt)
        snapdata = np.array([ t, r, rs, ld, sd, alpha, beta, gamma ])
        data.append(snapdata)

    data = np.array(data)
    df = pd.DataFrame({
        "Time / ns": data[:,0] / 1000.0,
        "r / A" : data[:,1],
        "rs / A" : data[:,2],
        "ld / A" : data[:,3],
        "sd / A" : data[:,4],
        "alpha / deg" : data[:,5],
        "beta / deg" : data[:,6],
        "gamma / deg" : data[:,7],
        })

    df.to_csv("%s.csv" % Opts["OutPre"], quoting=csv.QUOTE_NONNUMERIC, index=False)

    return df


if __name__ == '__main__':
    main()
