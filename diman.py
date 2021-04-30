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
    psi: float.
        Roll angle (in degrees).
    theta: float.
        Pitch angle (in degrees).
    phi: float.
        Yaw angle (in degrees).
    '''

    # Donor quantities
    # Invert order of principal axes so to have the long molecular axis first,
    # short second and the orthogonal one last, and check righthandedness.
    dcom = donor.center_of_mass()
    dpa = donor.principal_axes()[::-1].T
    if np.linalg.det(dpa) < 0:
        dpa[:,-1] = -dpa[:,-1]

    # Acceptor quantities
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
    psi, theta, phi = euler_angles_from_matrix(M)

    # Transform between -90 and 90
    psi = np.degrees(( np.radians(psi) + np.pi / 2 ) % np.pi - np.pi / 2)
    theta = np.degrees(( np.radians(theta) + np.pi / 2 ) % np.pi - np.pi / 2)
    phi = np.degrees(( np.radians(phi) + np.pi / 2 ) % np.pi - np.pi / 2)

    # Distance between coms
    r = acom - dcom
    rnorm = np.linalg.norm(r)

    # Project distance components onto the donor principal axes
    rdd, rsd, rpid = np.abs(np.dot(dpa, r))

    # Project distance components onto the accpt principal axes
    rda, rsa, rpia = np.abs(np.dot(apa, r))

    # Get the best of each pair
    rd = np.min([ rdd, rda ])
    rs = np.min([ rsd, rsa ])
    rpi = np.min([ rpid, rpia ])

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
        accpt = u.select_atoms(Opts["ASel"])
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
