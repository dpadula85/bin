#!/usr/bin/env python

import sys
import csv
import numpy as np
import pandas as pd
from tqdm import tqdm
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

    inp.add_argument('--top', type=str, dest='TopFile',
                     help='''Topology.''')

    inp.add_argument('--traj', type=str, dest='TrajFile',
                     help='''Trajectory.''')

    inp.add_argument('-s', '--sel', type=str, default=None, dest='Sel',
                     help='''Atom selection.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--outdir', default='.', type=str, dest='OutDir',
                     help='''Output File Format.''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


def compute_order_params(u):
    """
    Compute orientational order parameters for a set of molecular axes.

    Parameters
    ----------
    u : np.ndarray, shape (N, M, 3)
        Molecular principal axis of inertia. N is number of frames, M is number of residues.

    Returns
    -------
    P1 : np.ndarray, shape (N,)
        First-rank order parameter, averaged over residues.
    P2 : np.ndarray, shape (N,)
        Second-rank order parameter, derived from the Saupe tensor.
    P4 : np.ndarray, shape (N,)
        Fourth-rank order parameter, averaged over residues.
    alpha : np.ndarray, shape (N, M)
        Angles (in degrees) between molecular axis and director.
    n : np.ndarray, shape (N, 3)
        Director (principal orientation) vector for each frame.
    """
    # Compute Saupe tensor for each residue at each timestep
    Qs = 0.5 * (3 * np.einsum('ijk,ijl->ijkl', u, u) - np.eye(3)[None, None, :, :])

    # Average Saupe tensor over residues
    Q = Qs.mean(axis=1)

    # Diagonalise Saupe tensor to get eigenvalues and eigenvectors
    vals, vecs = np.linalg.eigh(Q)
    P2 = vals[:, -1]  # Largest eigenvalue corresponds to P2
    n = vecs[:, :, -1]  # Corresponding eigenvector is the director

    # Compute dot product between each molecular axis and the director
    dot = np.einsum('ijk,ijk->ij', u, n[:, None, :])
    norm_u = np.linalg.norm(u, axis=2)
    norm_n = np.linalg.norm(n, axis=1)
    dotp = dot / np.sqrt(norm_u * norm_n[:, None])

    # Compute angle alpha between axis and director
    alpha = np.arccos(np.clip(dotp, -1.0, 1.0))

    # First-rank order parameter: <cos(alpha)>
    P1 = np.cos(alpha).mean(axis=1)

    # Fourth-rank order parameter
    P4 = (0.125 * (35 * np.cos(alpha)**4 - 30 * np.cos(alpha)**2 + 3)).mean(axis=1)

    # Convert angle to degrees
    alpha = np.degrees(alpha)

    return P1, P2, P4, alpha, n


def main():
    """
    Main execution function for computing order parameters and saving results.

    Returns
    -------
    df : pd.DataFrame
        DataFrame containing P1, P2, P4 for full and core selections at each timestep.
    """
    args = options()
    u = mda.Universe(args['TopFile'], args['TrajFile'])

    if len(u.residues) == 0:
        raise ValueError("No residues found in the provided topology.")

    times = []
    n_frames = len(u.trajectory)
    n_residues = len(u.residues)

    us = np.zeros((n_frames, n_residues, 3))

    # Only allocate ucores if core selection was provided
    compute_core = args['Sel'] is not None
    ucores = np.zeros_like(us) if compute_core else None

    for i, ts in enumerate(tqdm(u.trajectory, desc="Processing frames")):
        times.append(ts.time)
        for j, res in enumerate(u.residues):
            pa = res.atoms.principal_axes()[::-1].T
            if np.linalg.det(pa) < 0:
                pa[:, -1] = -pa[:, -1]
            us[i, j] = pa[:, 0]

            if compute_core:
                core = res.atoms.select_atoms(args['Sel'])
                pacore = core.principal_axes()[::-1].T
                if np.linalg.det(pacore) < 0:
                    pacore[:, -1] = -pacore[:, -1]
                ucores[i, j] = pacore[:, 0]

    P1, P2, P4, alpha, n = compute_order_params(us)
    df_data = {
        "t": times,
        "P1": P1,
        "P2": P2,
        "P4": P4
    }

    df_alpha = pd.DataFrame(columns=["t"] + [f"resid_{i}" for i in range(1, n_residues + 1)],
                             data=np.c_[times, alpha])
    df_alpha.to_csv(f"{args['OutDir']}/alpha.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)

    if compute_core:
        P1c, P2c, P4c, alphac, nc = compute_order_params(ucores)
        df_data.update({
            "P1 core": P1c,
            "P2 core": P2c,
            "P4 core": P4c
        })
        df_alpha_core = pd.DataFrame(columns=["t"] + [f"resid_{i}" for i in range(1, n_residues + 1)],
                                     data=np.c_[times, alphac])
        df_alpha_core.to_csv(f"{args['OutDir']}/alpha_core.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)

    df = pd.DataFrame(df_data)
    df.to_csv(f"{args['OutDir']}/order.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)

    return df


if __name__ == "__main__":
    main()
