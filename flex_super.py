#!/usr/bin/env python

import numpy as np
import argparse as arg
from rdkit import Chem
from rdkit.Geometry import Point3D
from rdkit.Chem import AllChem, rdFMCS


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Docks a molecule to a reference.',
             formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")
    
    inp.add_argument('-r', '--ref', required=True, type=str, dest="RefSmiFile",
                     help='''Reference SMILES File.''')

    inp.add_argument('-c', '--coor', default=None, type=str, dest="RefXYZFile",
                     help='''Reference XYZ coordinate file.''')

    inp.add_argument('-n', '--new', required=True, type=str, dest="NewSmiFile",
                     help='''New Structure SMILES File.''')

    inp.add_argument('-m', '--map', default=None, type=str, dest="RefMapFile",
                     help='''Mapping between the reference SMILES and
                     coordinates.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--output', default="docked.xyz", type=str, dest='OutXYZFile',
                     help='''Output XYZ File Name.''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


def polish_smiles(smiles, kekule=False):
    '''
    Function to polish a SMILES string through the RDKit.

    Parameters
    ----------
    smiles: str.
        SMILES string.
    kekule: bool (default: False).
        whether to return Kekule SMILES.

    Returns
    -------
    polished: str.
        SMILES string.
    '''

    mol = Chem.MolFromSmiles(smiles)
    polished = Chem.MolToSmiles(mol, kekuleSmiles=kekule)

    return polished


def dock_to_ref(ref_smi, ref_xyz, new_smi, mapping=None):
    '''
    Function to dock a structure to a reference structure with coordinates.

    Parameters
    ----------
    ref_smi: str.
        Reference molecule SMILES string.
    ref_xyz: np.array (N, 3).
        Reference molecule coordinates.
    new_smi: str.
        New molecule (to be docked) SMILES string.
    mapping: dict (default: None).
        Dictionary mapping reference SMILES indices to reference coordinates.

    Returns
    -------
    coords: np.array (N, 4).
        XYZ structure of new molecule docked on top of the reference.
    '''

    # Create mol objs
    ref = Chem.MolFromSmiles(ref_smi)
    new = Chem.MolFromSmiles(new_smi)

    # Assign xyz coordinates to reference
    rwref = Chem.RWMol(ref)
    rwrefconf = Chem.Conformer(rwref.GetNumAtoms())
    for atidx in range(rwref.GetNumAtoms()):
        if mapping:
            xyz_idx = mapping[atidx]
        else:
            xyz_idx = atidx
        x, y, z = ref_xyz[xyz_idx]
        rwrefconf.SetAtomPosition(atidx, Point3D(x, y, z))

    # Find max common substruct mcs
    mcs = rdFMCS.FindMCS([ ref, new ])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    # Get portions corresponding to mcs
    ref_match = ref.GetSubstructMatch(mcs_mol)
    new_match = new.GetSubstructMatch(mcs_mol)

    # Assign coordinates to mcs
    rwmcs = Chem.RWMol(mcs_mol)
    rwconf = Chem.Conformer(rwmcs.GetNumAtoms())
    matches = rwmcs.GetSubstructMatch(mcs_mol)
    for i, atom in enumerate(matches):
        rwconf.SetAtomPosition(atom, rwrefconf.GetAtomPosition(ref_match[i]))

    rwmcs.AddConformer(rwconf)

    # Dock the new molecule to the mcs with the reference molecule
    AllChem.ConstrainedEmbed(new, rwmcs)

    # Get the final coordinates of the new molecule
    coords = []
    conf = new.GetConformer(0)
    for idx in range(new.GetNumAtoms()):
        z = new.GetAtomWithIdx(idx).GetAtomicNum()
        pos = conf.GetAtomPosition(idx)
        coord = [ z, pos.x, pos.y, pos.z ]
        coords.append(coord)

    coords = np.array(coords)

    return coords


def main(Opts):

    # Read reference
    with open(Opts["RefSmiFile"]) as f:
        smiles = f.readlines()[0].strip()
        ref_smi = polish_smiles(smiles)

    # Read coordinates and get rid of the atom column
    ref_xyz = np.genfromtxt(Opts["RefXYZFile"], skip_header=2)[:,1:]

    # Read mapping and make dict
    idxs = np.loadtxt(Opts["RefMapFile"]).astype(int)
    idxs_smi = (idxs[:,0] - 1).tolist()
    idxs_xyz = (idxs[:,1] - 1).tolist()
    smi_to_xyz = dict(zip(idxs_smi, idxs_xyz))

    # Read new molecule
    with open(Opts["NewSmiFile"]) as f:
        smiles = f.readlines()[0].strip()
        new_smi = polish_smiles(smiles)

    # Dock
    new_xyz = dock_to_ref(ref_smi, ref_xyz, new_smi, mapping=smi_to_xyz)

    # Write result to output
    with open(Opts["OutXYZFile"], "w") as f:
        f.write("%d\n\n" % len(new_xyz))
        np.savetxt(f, new_xyz, fmt="%-5d %12.6f %12.6f %12.6f")

    return


if __name__ == '__main__':
    Opts = options()
    main(Opts)
