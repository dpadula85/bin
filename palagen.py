#!/usr/bin/env python

# Splits a protein in residues for EET calculations.
# The ith residue is described as the N-MethylAmide of Alanine, as follows:
# it includes all the atoms of the ith residue and the NH and CA, HA of the
# ith + 1 residue. The N of the ith residue is capped by keeping the carbonyl
# C of the ith - 1 residue, and transforming it into an H at 1.08 A distance.
# To split the protein in each point it invokes qmip.py by Lorenzo, which
# splits a structure for a QM/MMPol calculation. In this case, we will ignore
# the MMPol part, and keep only the QM one. The procedure is repeated for
# each residue in the protein.

# Import standard Python modules
import os
import subprocess
import shutil
import numpy as np

# Import personal modules
import G09_files

# Read mol2 file adapted from Lorenzo's qmip
def parse_mol2(mol2file):

    with open(mol2file) as f:
        FoundAt = False
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
                try:
                    NRes = int(info[2])
                except:
                    NRes = 1

                # Lorenzo's code
                #AtNames = []
                #AtTypes = []
                #AtRes   = []
                #AtResId = []
                #AtCoord = []
           
                # Same variables as Lorenzo, I just don't like his names
                atom_names = []
                atom_types = []
                res_names = []
                res_ids = []
                atom_coord = []

            # Read Atoms
            elif  line[0:13] == '@<TRIPOS>ATOM':
                for i in range(NAtoms):
                    data = f.readline().split()
                    
                    # Lorenzo's code
                    #AtNames += [data[1]]
                    #AtTypes += [data[5]]
                    #AtRes += [data[7]]
                    #AtCoord += data[2:5]
                    #RID = int(data[6])
                    #AtResId += [RID]

                    # I don't like Lorenzo's structure of the data
                    # since it is just a plain list of atoms, residues
                    # and ids. Residues and ids are repeated for each atom
                    # in the residue, thus it is very redundant.
                    # The following code creates a different structure:
                    # res_names and res_ids are plain lists, and the elements
                    # appear only once per residue.
                    # atom_names, type and coord are lists of lists. Each sublist
                    # corresponds to a residue and contains all the info
                    # for  on the atoms in that residue.
                    # This way it is also easy to verify that the first five
                    # atoms in each residue are N, CA, C, O, CB.

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
                        atom_coord.append(coords)
                    
                        # initialize data for the new one
                        names = [data[1]]
                        types = [data[5]]
                        coords = [data[2:5]]
                    
                    # if still in the old residue
                    else:
                    
                        # save data for the new atom
                        try:
                            names.append(data[1])
                            types.append(data[5])
                            coords.append(data[2:5])
                    
                        # unless no atom has been saved before
                        except:
                            names = [data[1]]
                            types = [data[5]]
                            coords = [data[2:5]]

                FoundAt = True
                if FoundAt:
                    # save data for the last residue 
                    atom_names.append(names)
                    atom_types.append(types)
                    atom_coord.append(coords)


    return atom_names, atom_types, res_names, res_ids
    #return AtNames, AtTypes, AtRes, AtResId

# Now I gotta use the information parsed from the mol2 file to set up the
# selection.in file for a correct cut in qmip.
# Keep in mind that N, CA, C, O, CB should always be the first 5 atoms and that
# the first four H atoms in a residue are the ones connected respectively
# to N (one H), CA (one H), CB (two Hs). These are the atoms I wanna use
# to set up a split and link in qmip.

#
# =========================
#  The Program Starts Here
# =========================
#

if __name__ == '__main__':
    
    # qmip general options
    qmip_db_file ='full.db' 
    qmip_mol2_file = 'prova.mol2'
    qmip_infile = 'selection.in'
    qmip_cmd = 'qmip.py --db %s --mol2 %s -i %s' % (qmip_db_file, qmip_mol2_file, qmip_infile)
    qmip_outfile = 'frame00001.com'

    # parse mol2 files
    atom_names, atom_types, res_names, res_ids = parse_mol2(qmip_mol2_file)

    # associate res_names to res_ids
    prot_seq = dict(zip(res_ids, res_names))

    # associate each aminoacid to a list of atoms 
    aa_atoms = dict(zip(res_names, atom_names))

    # associate each aminoacid to some cut options 
    split_dict = {}
    split_dict['ALA'] = '' 
    split_dict['ARG'] = '1-5,12-15 link 5 6\n'
    split_dict['ASN'] = '1-5,9-12 link 5 6\n'
    split_dict['ASP'] = '1-5,9-12 link 5 6\n'
    split_dict['CYS'] = '1-5,7-10 link 5 6\n'
    split_dict['GLN'] = '1-5,10-13 link 5 6\n'
    split_dict['GLU'] = '1-5,10-13 link 5 6\n'
    split_dict['GLY'] = ''
    split_dict['HIS'] = '1-5,11-14 link 5 6\n'
    split_dict['ILE'] = '1-5,9-11 link 5 6 5 7\n'
    split_dict['LEU'] = '1-5,9-12 link 5 6\n'
    split_dict['LYS'] = '1-5,10-13 link 5 6\n'
    split_dict['MET'] = '1-5,9-12 link 5 6\n'
    split_dict['PHE'] = '1-5,12-15 link 5 6\n'
    split_dict['PRO'] = '1-5,8-10 link 1 7 5 6\n'
    split_dict['SER'] = '1-5,7-10 link 5 6\n'
    split_dict['THR'] = '1-5,8-10 link 5 6 5 7\n'
    split_dict['TRP'] = '1-5,15-18 link 5 6\n'
    split_dict['TYR'] = '1-5,13-16 link 5 6\n'
    split_dict['VAL'] = '1-5,8-10 link 5 6 5 7\n'

    # General Options to set up G09 input file for eet=prop calculation on the residue
    opts = {}
    opts['nproc'] = 12
    opts['mem'] = 2500
    opts['mem_unit'] = 'MW'
    opts['funct'] = 'cam-b3lyp'
    opts['basis'] = '6-31G(d)'
    opts['job'] = 'td=nstates=10 eet=prop IOp(2/12=3) scrf=(iefpcm,solvent=water)'
    opts['add_opts'] = '\ng03defaults nocav nodis norep rmin=0.5 ofac=0.8 sphereonh=4\n\n' 

    # write coordinates for each residue in an .xyz file
    xyz = open('fragment.xyz', 'w')

    for res_id, res_name in prot_seq.iteritems():
        
        # special case for the last aminoacid
        # which should be skipped
        if res_id == max(prot_seq.keys()):
            continue

        # special case for the first aminoacid
        # the previous residue should not be included
        if res_id == min(prot_seq.keys()):

            # qmip residue-specific options
            with open(qmip_infile, 'w') as sel_f:
            
                # two residues to include
                sel_f.write('\nqmresid %d %d\n' % (res_id, res_id + 1))

                # splitting fot the ith residue 
                sel_f.write('split %d  1-5,10-14 link 5 6\n' % res_id )
            
                # for the ith + 1 residue keep always N, CA (atoms 1, 2) and
                # the first two H atoms (H and HA). BEWARE!!! GLY and PRO are
                # special cases, so I gotta handle them separately.

                if prot_seq[res_id + 1] == 'GLY':
                    sel_f.write('split %d 1,2,5-7 link 2 3\n\n' % (res_id + 1))

                elif prot_seq[res_id + 1] == 'PRO':
                    sel_f.write('split %d 1,2,8 link 1 7 2 3\n\n' % (res_id + 1))

                else:
                    # Get the indexes for the first two H atoms
                    target = 2
                    count = 0
                    H_idxs = []
                    for idx, atom in enumerate(aa_atoms[prot_seq[res_id + 1]], start=1):
                        if atom == 'H' or atom =='HA':
                            H_idxs.append(idx)
                            count += 1
                        if count == target:
                            break
                    
                    # assign the indexes to some variable for easier writing
                    H, HA = H_idxs
                    sel_f.write('split %d 1,2,%d,%d link 2 3 2 5\n\n' % (res_id + 1, H, HA))

        # in any other case
        else:

            # qmip residue-specific options
            with open(qmip_infile, 'w') as sel_f:
            
                # three residues to include
                sel_f.write('\nqmresid %d %d %d\n' % (res_id - 1, res_id, res_id + 1))
            
                # keep carbonylic C of the ith - 1 residue, which is always atom 3
                sel_f.write('\nsplit %d 3\n' % (res_id - 1))
            
                # write splitting options for the current residue according to the
                # splitting dictionary
                if res_name == 'ALA' or res_name == 'GLY':
                    pass

                else:
                    sel_f.write('split %d %s' % (res_id, split_dict[res_name]))
            
                # for the ith + 1 residue keep always N, CA (atoms 1, 2) and
                # the first two H atoms (H and HA). BEWARE!!! GLY and PRO are
                # special cases, so I gotta handle them separately.

                if prot_seq[res_id + 1] == 'GLY':
                    sel_f.write('split %d 1,2,5-7 link 2 3\n' % (res_id + 1))

                elif prot_seq[res_id + 1] == 'PRO':
                    sel_f.write('split %d 1,2,8 link 1 7 2 3\n\n' % (res_id + 1))

                else:
                    # Get the indexes for the first two H atoms
                    target = 2
                    count = 0
                    H_idxs = []
                    for idx, atom in enumerate(aa_atoms[prot_seq[res_id + 1]], start=1):
                        if atom == 'H' or atom =='HA':
                            H_idxs.append(idx)
                            count += 1
                        if count == target:
                            break
                    
                    # assign the indexes to some variable for easier writing
                    H, HA = H_idxs
                    sel_f.write('split %d 1,2,%d,%d link 2 3 2 5\n\n' % (res_id + 1, H, HA))
        
        # execute qmip
        subprocess.call(qmip_cmd.split())
        #print res_id, res_name
        
        # after qmip executes, a file called frame00001.com is generated for our
        # further processing. Let's open it with my G09 module.
        bad_f = G09_files.input_file(qmip_outfile)
        
        # First of all, the C carbon of residue i - 1 shall be substituted by an
        # H atom. To do this, we need the coordinates of this C atom and of the N
        # atom to be capped with an H. These two atoms are the first two, since
        # C in is in residue i - 1 and N is the first atom of the ith residue
        # since the protein starts from the N-terminal side
        if res_id > min(prot_seq.keys()):
 
            bad_C = np.array(bad_f.structure[0][1:])
            cap_N = np.array(bad_f.structure[1][1:])
            
            # To substitute the C with an H atom, let's define the versor of the
            # C-N bond, and multiply it for the typical N-H distance in A.
            d_CN = bad_C - cap_N
            v = d_CN / np.linalg.norm(d_CN)
            d_NH = v * 0.86 + cap_N
            cap_H = ['1'] + d_NH.tolist()
            
            # substitute the C atom with the just generated H one
            bad_f.structure[0] = cap_H
        
        # Set residue-specific options for eet=prop calculation
        opts['chk'] = '%04d' % (res_id + 1)
        opts['title'] = '%04d properties calculation' % (res_id + 1)
        opts['charge'] = bad_f.charge
        opts['mult'] = bad_f.mult
        opts['structure'] = bad_f.structure
        
        # write the file with the correctly capped structure and options
        prop_infile = G09_files.input_file('%04d.com' % (res_id + 1), opts)

        # move the just written file to its directory
        os.makedirs(os.path.join(os.getcwd(), '%04d' % (res_id + 1)))
        shutil.move(prop_infile.name, os.path.join(os.getcwd(), '%04d' % (res_id + 1)))

        # write coordinates for this residue in the xyz file
        for xyz_coor in prop_infile.structure:
            xyz.write('%s %12.8f %12.8f %12.8f\n' % tuple(xyz_coor))
        pass

    xyz.close()
