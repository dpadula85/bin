#!/usr/bin/env python

# Import standard Python modules
import re
import itertools
import os
import shutil
import numpy as np
from difflib import SequenceMatcher

# Import personal modules
import DFT_methods as dft
from elements import ELEMENTS
import util as u

#
# ================
#  G09 input file
# ================
#


class input_file:
    '''A class describing G09 input files.'''

    def __init__(self, infile, opts_dict=None):

        self.name = os.path.split(infile)[1]
        self.path = os.path.split(infile)[0]
        self.file = os.path.join(self.path, self.name)

        if not os.path.exists(infile) or os.stat(infile).st_size == 0:
            self.options = self.set_options(opts_dict)
            self.write_job()

        if os.path.exists(infile) and os.stat(infile).st_size > 0:
            self.options = self.get_options()

    def set_options(self, options=None):
        '''Default options for creation of a new file.'''

        default = {'name': self.name,
                   'nproc': 16,
                   'mem': 24,
                   'mem_unit': 'GB',
                   'chk': self.name,
                   'funct': 'b3lyp',
                   'basis': '6-31g*',
                   'job': 'td=(nstates=16)',
                   'title': '%s td' % self.name,
                   'charge': 0,
                   'mult': 1,
                   'structure': [['O', 0.000000, 0.000000, 0.000000],
                                 ['H', 0.758602, 0.000000, 0.504284],
                                 ['H', 0.758602, 0.000000, -0.504284]],
                   'add_opts': '\n'}

        if not options:
            options = default
        else:
            options = u.dict_compare(default, options)

        # Needed in case another option dictionary is passed to the function
        options['name'] = self.name

        return options

    def get_options(self):

        self.nproc = self.get_nproc()
        self.mem, self.mem_unit = self.get_mem()
        self.chk = self.get_chk()
        self.keywords = self.get_keywords()
        self.funct = self.get_funct()
        self.basis = self.get_basis()
        self.job = self.get_job()
        self.title = self.get_title()
        self.charge = self.get_charge()
        self.mult = self.get_mult()
        self.structure = self.get_structure()
        self.add_opts = self.get_add_opts()
        self.natoms = len(self.structure)

        options = {'name': self.name,
                   'nproc': self.nproc,
                   'mem': self.mem,
                   'mem_unit': self.mem_unit,
                   'chk': self.chk,
                   'funct': self.funct,
                   'basis': self.basis,
                   'job': self.job,
                   'title': self.title,
                   'charge': self.charge,
                   'mult': self.mult,
                   'structure': self.structure,
                   'add_opts': self.add_opts}

        return options

    def __str__(self):
        '''Prints a string containing the whole file.'''

        with open(self.file, 'r') as f:
            print(f.read())

    def get_nproc(self):
        '''Returns the number of processors.'''

        with open(self.file, 'r') as f:
            for line in f:
                if 'nproc' in line:
                    nproc = int(line.split('=')[-1].strip())
                    break
                else:
                    nproc = None

        return nproc

    def get_mem(self):
        '''Returns the memory and its unit.'''

        with open(self.file, 'r') as f:
            for line in f:
                if 'mem' in line:
                    mem_str = line.split('=')[-1].strip()
                    mem, unit = re.search('(\d+)([a-zA-Z]+)', mem_str).groups()
                    break
                else:
                    mem, unit = None, ''

        if mem:
            mem = int(mem)

        return mem, unit

    def get_chk(self):
        '''Returns the name of the checkpoint file without extension.'''

        with open(self.file, 'r') as f:
            for line in f:
                if 'chk' in line:
                    chk_str = line.split('=')[-1].strip()
                    chk = chk_str.split('.')[0]
                    break
                else:
                    chk = ''

        return chk

    def get_keywords(self):
        '''Returns the full string of keywords.'''

        with open(self.file, 'r') as f:
            for line in f:
                if '#' in line:
                    keywords = line.strip()
                    break
                else:
                    keywords = ''

        return keywords

    def get_funct(self):
        '''Returns the functional.'''

        # The order of the keywords is not fixed. This function parses the
        # keywords and compares it to the available functionals in G09.

        if self.keywords:
            for keyword in self.keywords.split():
                if keyword.lower() in map(
                        str.lower, list(itertools.chain.from_iterable(
                            dft.functionals.values()))):
                    funct = keyword
                    break
                else:
                    funct = ''

        return funct

    def get_basis(self):
        '''Returns the basis set.'''

        # The order of the keywords is not fixed. This function parses the
        # keywords and compares it to the available basis sets in G09. Since
        # slight modifications in the names can occur due to diffuse or
        # polarized functions, a 65% similarity threshold is used to determine
        # whether a keyword is a basis set.

        stop = False
        if self.keywords:
            for keyword in self.keywords.split():
                for basis_set in map(
                        str.lower, list(itertools.chain.from_iterable(
                            dft.basis.values()))):
                    if SequenceMatcher(
                            None, keyword.lower(), basis_set).ratio() > 0.65:
                        basis = keyword
                        stop = True
                        break
                if stop:
                    break
                else:
                    basis = ''

        return basis

    def get_job(self):
        '''Returns the job-related keywords.'''

        # Once the functional and the basis set are determined, the job-related
        # keywords are obtained as keywords except functional and basis set.

        if self.keywords:

            # Starting from 1 to get rid of the starting '#' or '#p'
            keywords = self.keywords.split()[1:]
            job = [x for x in keywords if x != self.funct and x != self.basis]

        return ' '.join(job)

    def get_title(self):
        '''Returns the title.'''

        with open(self.file, 'r') as f:
            for line in f:
                if not line.strip():
                    title = next(f).strip()
                    break
                else:
                    title = ''

        return title

    def get_charge(self):
        '''Returns the charge.'''

        # Pattern describing the charge-multiplicity line
        pattern = re.compile('(-?\d\s\d)')

        with open(self.file, 'r') as f:
            for line in f:
                if pattern.match(line):
                    charge = line.split()[0]
                    break
                else:
                    charge = None

        if charge:
            charge = int(charge)

        return charge

    def get_mult(self):
        '''Returns the multiplicity.'''

        # Pattern describing the charge-multiplicity line
        # pattern = re.compile('(\d\s\d)')
        pattern = re.compile('(-?\d\s\d)')

        with open(self.file, 'r') as f:
            for line in f:
                if pattern.match(line):
                    mult = line.split()[1]
                    break
                else:
                    mult = None

        if mult:
            mult = int(mult)

        return mult

    def get_structure(self):
        '''Returns the structure.'''

        structure = []
        opt = 0

        # Pattern describing the charge-multiplicity line
        # pattern = re.compile('(\d\s\d)')
        pattern = re.compile('(-?\d\s\d)')

        with open(self.file, 'r') as f:
            for line in f:
                if pattern.match(line):
                    opt = 1

                if not line.strip():
                    opt = 0

                if opt == 1 and not pattern.match(line):
                    curr_line = filter(None, line.split())
                    atom = curr_line[0]
                    atom_x = float(curr_line[1])
                    atom_y = float(curr_line[2])
                    atom_z = float(curr_line[3])
                    structure.append([atom, atom_x, atom_y, atom_z])

        return structure

    def get_add_opts(self):
        '''Returns additional options.'''

        # Usually additional options are just a blank line. In many cases
        # it's useful to retrieve it (additional functions in the basis set,
        # solvation spheres etc).

        # Pattern describing the charge-multiplicity line
        # pattern = re.compile('(\d\s\d)')
        pattern = re.compile('(-?\d\s\d)')

        opt = 0
        switch = 0
        add_opts = []

        with open(self.file, 'r') as f:
            for line in f:
                if pattern.match(line):
                    opt = 1

                if opt == 1:
                    switch = 1

                if not line.strip():
                    opt = 0

                if opt == 0 and switch == 1:
                    add_opts.append(line)

        return add_opts

    def get_com(self):
        '''Returns the coordinates of the center of mass.'''

        tot_mass = 0
        x_com = 0
        y_com = 0
        z_com = 0

        for atom in self.structure:

            try:
                m = ELEMENTS[atom[0]].mass
            except:
                m = ELEMENTS[int(atom[0])].mass

            tot_mass += m
            x_com += atom[1] * m
            y_com += atom[2] * m
            z_com += atom[3] * m

        x_com = x_com / tot_mass
        y_com = y_com / tot_mass
        z_com = z_com / tot_mass

        return (x_com, y_com, z_com)

    def write_job(self, opts_dict=None):
        '''Writes a G09 input for the next calculation. Options for the new
        job can be specified, otherwise default options will be used.'''

        if not opts_dict:
            opts_dict = self.options

        if opts_dict:
            opts_dict = u.dict_compare(self.options, opts_dict)

        with open(opts_dict['name'], 'w') as f:

            if opts_dict['nproc']:
                f.write('%%nproc=%d\n' % opts_dict['nproc'])

            if opts_dict['mem']:
                f.write('%%mem=%d%s\n' %
                        (opts_dict['mem'], opts_dict['mem_unit']))

            if opts_dict['chk']:
                f.write('%%chk=%s.chk\n' % opts_dict['chk'])

            f.write('#p %s %s %s\n' %
                    (opts_dict['funct'], opts_dict['basis'], opts_dict['job']))

            f.write('\n')
            f.write('%s\n' % opts_dict['title'])
            f.write('\n')
            f.write('%d %d\n' % (opts_dict['charge'], opts_dict['mult']))

            for atom in opts_dict['structure']:
                f.write('%2s\t%12.8f\t%12.8f\t%12.8f\n' %
                        (atom[0], atom[1], atom[2], atom[3]))

            for opt in opts_dict['add_opts']:
                f.write(opt)

#
# ================
#  G09 output file
# ================
#
# Cannot inherit from G09 input class since methods are slightly different
# due to small differences in the information format.
# Methods to retrieve basic information about the file are the same, corrected
# because of the format differences.
# Additional methods related to the job performed are written in child classes
# to be defined afterwards
#


class output_file:
    '''A class describing G09 logfiles.'''

    def __init__(self, outfile):

        self.name = os.path.split(outfile)[1]
        self.path = os.path.split(outfile)[0]
        self.file = os.path.join(self.path, self.name)
        self.options = self.get_options()

    def get_options(self):

        self.nproc = self.get_nproc()
        self.mem, self.mem_unit = self.get_mem()
        self.chk = self.get_chk()
        self.keywords = self.get_keywords()
        self.funct = self.get_funct()
        self.basis = self.get_basis()
        self.job = self.get_job()
        self.title = self.get_title()
        self.charge = self.get_charge()
        self.mult = self.get_mult()
        self.structure = self.get_initial_structure()
        self.add_opts = '\n'
        self.natoms = len(self.structure)
        self.types = self.get_types()

        options = {'name': self.name,
                   'nproc': self.nproc,
                   'mem': self.mem,
                   'mem_unit': self.mem_unit,
                   'chk': self.chk,
                   'funct': self.funct,
                   'basis': self.basis,
                   'job': self.job,
                   'title': self.title,
                   'charge': self.charge,
                   'mult': self.mult,
                   'structure': self.structure,
                   'add_opts': self.add_opts}

        return options

    def __str__(self):
        '''Prints a string containing the whole file.'''

        with open(self.file, 'r') as f:
            print(f.read())

    def get_types(self):
        '''Sets a label for the performed G09 job.'''

        # This function is actually not used yet.
        # The idea is to determine what kind of job was run
        # and create an instance of a child class based on it.

        keywords_split = re.split('[ =()]', self.job)
        job_types = [dft.job_types[keyword]
                     for keyword in keywords_split
                     if keyword in dft.job_types.keys()]

        return job_types

    def get_end(self):
        '''Checks if the G09 job terminated correctly.'''

        with open(self.file, 'r') as f:
            for line in f:
                if "Normal termination" in line:
                    status = 'Finished'
                    break
                else:
                    status = 'Error'

        return status

    def get_nproc(self):
        '''Returns the number of processors.'''

        with open(self.file, 'r') as f:
            for line in f:
                if 'nproc' in line.lower():
                    nproc = int(line.split('=')[-1].strip())
                    break
                else:
                    nproc = None

        return nproc

    def get_mem(self):
        '''Returns the memory and its unit.'''

        with open(self.file, 'r') as f:
            for line in f:
                if '%mem' in line.lower():
                    mem_str = line.split('=')[-1].strip()
                    mem, unit = re.search('(\d+)([a-zA-Z]+)', mem_str).groups()
                    break
                else:
                    mem, unit = None, ''

        if mem:
            mem = int(mem)

        return mem, unit

    def get_chk(self):
        '''Returns the name of the checkpoint file without extension.'''

        with open(self.file, 'r') as f:
            for line in f:
                if 'chk' in line.lower():
                    chk_str = line.split('=')[-1].strip()
                    chk = chk_str.split('.')[0]
                    break
                else:
                    chk = ''

        return chk

    def get_keywords(self):
        '''Returns the full string of keywords.'''

        # Annoying syntax with switch because G09 splits the keywords
        # on more lines if they are too many

        previous = None
        switch = 0
        temp = []

        with open(self.file, 'r') as f:
            for line in f:
                if '#' in line and '----' in previous:
                    switch = 1

                if '----' in line:
                    switch = 0

                if switch == 1:
                    temp.append(line.strip())

                previous = line

        if temp:
            keywords = ''.join(temp)
        else:
            keywords = ''

        return keywords

    def get_funct(self):
        '''Returns the functional.'''

        # The order of the keywords is not fixed. This function parses the
        # keywords and compares it to the available functionals in G09.

        if self.keywords:
            for keyword in self.keywords.split():
                if keyword.lower() in map(
                        str.lower, list(itertools.chain.from_iterable(
                            dft.functionals.values()))):
                    funct = keyword
                    break
                else:
                    funct = ''

        return funct

    def get_basis(self):
        '''Returns the basis set.'''

        # The order of the keywords is not fixed. This function parses the
        # keywords and compares it to the available basis sets in G09. Since
        # slight modifications in the names can occur due to diffuse or
        # polarized functions, a 65% similarity threshold is used to determine
        # whether a keyword is a basis set.

        stop = False
        if self.keywords:
            for keyword in self.keywords.split():
                for basis_set in map(
                        str.lower, list(itertools.chain.from_iterable(
                            dft.basis.values()))):
                    if SequenceMatcher(
                            None, keyword.lower(), basis_set).ratio() > 0.65:
                        basis = keyword
                        stop = True
                        break
                if stop:
                    break
                else:
                    basis = ''

        return basis

    def get_job(self):
        '''Returns the job-related keywords.'''

        # Once the functional and the basis set are determined, the job-related
        # keywords are obtained as keywords except functional and basis set.

        if self.keywords:
            keywords = self.keywords.split()[1:]
            job = [x for x in keywords if x != self.funct and x != self.basis]

        return ' '.join(job)

    def get_title(self):
        '''Returns the title.'''

        preprevious = None
        previous = None
        with open(self.file, 'r') as f:
            for line in f:
                if "Symbolic Z-matrix" in line:
                    title = preprevious.strip()
                    break
                else:
                    title = ''
                preprevious, previous = previous, line

        return title

    def get_charge(self):
        '''Returns the charge.'''

        with open(self.file, 'r') as f:
            for line in f:
                if "multiplicity" in line.lower():
                    charge = line.split()[2]
                    break
                else:
                    charge = None

        if charge:
            charge = int(charge)

        return charge

    def get_mult(self):
        '''Returns the multiplicity.'''

        with open(self.file, 'r') as f:
            for line in f:
                if "multiplicity" in line.lower():
                    mult = line.split()[5]
                    break
                else:
                    mult = None

        if mult:
            mult = int(mult)

        return mult

    def get_initial_structure(self):
        '''Returns the structure written in the corresponding input.'''

        # Adding previous line control to avoid crash when multiple jobs are
        # executed in chain
        previous = None
        structure = []
        opt = 0

        with open(self.file, 'r') as f:
            for line in f:
                if "multiplicity" in line.lower() and "Symbolic Z-matrix" in \
                        previous:
                    opt = 1

                if not line.strip():
                    opt = 0

                if opt == 1 and "multiplicity" not in line.lower():
                    curr_line = filter(None, line.split())
                    atom = curr_line[0]
                    atom_x = float(curr_line[1])
                    atom_y = float(curr_line[2])
                    atom_z = float(curr_line[3])
                    structure.append([atom, atom_x, atom_y, atom_z])

                previous = line

        return structure

    def get_com(self):
        '''Returns the coordinates of the center of mass.'''

        tot_mass = 0
        x_com = 0
        y_com = 0
        z_com = 0

        for atom in self.structure:

            try:
                m = ELEMENTS[atom[0]].mass
            except:
                m = ELEMENTS[int(atom[0])].mass

            tot_mass += m
            x_com += atom[1] * m
            y_com += atom[2] * m
            z_com += atom[3] * m

        x_com = x_com / tot_mass
        y_com = y_com / tot_mass
        z_com = z_com / tot_mass

        return (x_com, y_com, z_com)

    def set_options(self, job=None):
        '''Sets the default options for a G09 job.'''

        if not self.add_opts:
            self.add_opts = '\n'

        if not job and 'opt' in self.types:
            self.job = 'td=(nstates=16) nosymm'

        options = {'name': self.name,
                   'nproc': self.nproc,
                   'mem': self.mem,
                   'mem_unit': self.mem_unit,
                   'chk': self.chk,
                   'funct': self.funct,
                   'basis': self.basis,
                   'job': self.job,
                   'title': self.title,
                   'charge': self.charge,
                   'mult': self.mult,
                   'structure': self.structure,
                   'add_opts': self.add_opts}

        return options

    def write_next_job(self, opts_dict=None):
        '''Writes a G09 input for the next calculation. Options for the new
        job can be specified, otherwise default options will be used.'''

        if not opts_dict:
            opts_dict = self.options

        if opts_dict:
            opts_dict = u.dict_compare(self.options, opts_dict)

        input_file('%s_new.com' % opts_dict['name'].split('.')[0], opts_dict)

#
# ==============================
#  G09 optimization output file
# ==============================
#
# Inherits from G09 output class. This class should contain methods specific
# to extract data from G09 optimizations.
#


class opt_output_file(output_file):
    '''A class describing G09 optimization logfiles.'''

    def __init__(self, outfile):
        output_file.__init__(self, outfile)
        self.results = self.get_results()

    def get_opt_structure(self):
        '''Returns the optimized structure.'''

# Example:
#
#                         Standard orientation:
# ---------------------------------------------------------------------
# Center     Atomic      Atomic             Coordinates (Angstroms)
# Number     Number       Type             X           Y           Z
# ---------------------------------------------------------------------
#      1          7           0       -0.904158   -2.566751    0.918789

        structure = []
        opt = 0
        with open(self.file, 'r') as f:
            for line in f:
                if "Optimization completed" in line:
                    opt = 1

                # It could be Standard or Input orientation
                # depending on nosymm keyword
                if "orientation:" in line and opt == 1:
                    # Skip the head-of-table lines (4 lines)
                    next(f)
                    next(f)
                    next(f)
                    next(f)
                    for i in range(self.natoms):
                        curr_line = next(f).split()
                        atom_n = int(curr_line[1])
                        atom = ELEMENTS[atom_n].symbol
                        atom_x = float(curr_line[3])
                        atom_y = float(curr_line[4])
                        atom_z = float(curr_line[5])
                        structure.append([atom, atom_x, atom_y, atom_z])
                    opt = 0

        return structure

    def get_opt_en(self):
        '''Returns the energy of the optimized structure.'''

        with open(self.file, 'r') as f:
            for line in f:
                if "SCF Done" in line:
                    energy = float(line.split()[4])

        return energy

    def get_results(self):
        '''Returns the results of the optimization.'''

        self.structure = self.get_opt_structure()
        self.energy = self.get_opt_en()
        self.energy_unit = 'au'

        results = {'optimized structure': self.structure,
                   'optimized energy': self.energy,
                   'energy unit': self.energy_unit}

        return results

    def en_convert(self, new_u):
        '''Converts the energy of the optimized structure to another unit.'''

        u1 = self.energy_unit
        u2 = new_u
        factor1 = u.energy_conversion[u1]
        factor2 = u.energy_conversion[u2]

        new_energy = self.energy * factor2 / factor1

        self.energy_unit = u2
        self.energy = new_energy

        print('%f %s' % (self.energy, self.energy_unit))

#
# ====================
#  G09 td output file
# ====================
#
# Inherits from G09 output class. This class should contain methods specific
# to extract data from G09 td calculations.
#


class td_output_file(output_file):
    '''A class describing G09 td logfiles.'''

    def __init__(self, outfile):
        output_file.__init__(self, outfile)
        self.results = self.get_results()

    def get_es_prop(self):
        '''Extracts excitation energies and wavelengths.'''

        energies = []
        wavelengths = []
        oscillators = []

        with open(self.file, 'r') as f:
            for line in f:

                if "Excited State  " in line:

# Example:
#
# Excited State   1:      Singlet-?Sym    4.8743 eV  254.36 nm  f=0.0151  <S**2>=0.000

                    energy = line.split()[4]
                    energies.append(float(energy))

                    wavelength = line.split()[6]
                    wavelengths.append(float(wavelength))

                    f_osc = line.split()[-2].split('=')[-1]
                    oscillators.append(float(f_osc))

        return energies, wavelengths, oscillators

    def get_es_rot(self):
        '''Extracts rotatory strenghts.'''

        rot_vel = []
        rot_len = []

        with open(self.file, 'r') as f:
            for line in f:

                if "R(velocity)" in line:

# Example:
#
#      state          XX          YY          ZZ    R(velocity)    E-M Angle
#        1       -143.5843   -117.0927    -17.9532    -92.8767      136.17

                    for i in range(self.nstates):
                        r_vel = next(f).split()[-2]
                        rot_vel.append(float(r_vel))

                if "R(length)" in line:

# Example:
#
#       state          XX          YY          ZZ       R(length)
#         1        -133.6881     23.8640   -170.4743    -93.4328

                    for i in range(self.nstates):
                        r_len = next(f).split()[-1]
                        rot_len.append(float(r_len))

        return rot_vel, rot_len

    def get_mu(self):
        '''Extracts transition electric dipoles.'''

        mu_vel = []
        mu_len = []

        with open(self.file, 'r') as f:
            for line in f:

                if "excited state transition electric dipole" in line:

# Example:
#
# Ground to excited state transition electric dipole moments (Au):
#       state          X           Y           Z        Dip. S.      Osc.
#         1        -0.8862      0.3776     -0.1925      0.9651      0.1142

                    next(f)
                    for i in range(self.nstates):
                        curr_line = next(f)
                        mu_len_x = float(curr_line.split()[1])
                        mu_len_y = float(curr_line.split()[2])
                        mu_len_z = float(curr_line.split()[3])
                        mu_len.append([mu_len_x, mu_len_y, mu_len_z])

                if "excited state transition velocity dipole" in line:

# Example:
#
# Ground to excited state transition velocity dipole moments (Au):
#       state          X           Y           Z        Dip. S.      Osc.
#         1         0.1548     -0.0673      0.0468      0.0307      0.1152

                    next(f)
                    for i in range(self.nstates):
                        curr_line = next(f)
                        mu_vel_x = float(curr_line.split()[1])
                        mu_vel_y = float(curr_line.split()[2])
                        mu_vel_z = float(curr_line.split()[3])
                        mu_vel.append([mu_vel_x, mu_vel_y, mu_vel_z])

        return mu_vel, mu_len

    def get_mag(self):
        '''Extracts transition magnetic dipoles.'''

        m = []

        with open(self.file, 'r') as f:
            for line in f:

                if "excited state transition magnetic dipole" in line:

# Example:
#
# Ground to excited state transition magnetic dipole moments (Au):
#       state          X           Y           Z
#         1         5.9417      4.8361    -12.7233

                    next(f)
                    for i in range(self.nstates):
                        curr_line = next(f)
                        m_x = float(curr_line.split()[1])
                        m_y = float(curr_line.split()[2])
                        m_z = float(curr_line.split()[3])
                        m.append([m_x, m_y, m_z])

        return m

    def mag_int(self, s):
        '''Returns the intrinsic component of the transition magnetic dipole
        moment relative to state s.'''

        # m = m_int + r x mu
        x, y, z = self.get_com()
        r = np.array([x, y, z])
        m = np.array(self.m[s - 1])
        mu = np.array(self.mu_len[s - 1])
        
        mag_int = m - np.cross(r, mu)

        return mag_int.tolist()

    def get_results(self):
        '''Returns the results of the excited states calculation.'''

        self.es_energies, self. es_wavelengths, self.f_osc = self.get_es_prop()
        self.nstates = len(self.es_energies)
        self.r_vel, self.r_len = self.get_es_rot()
        self.mu_vel, self.mu_len = self.get_mu()
        self.m = self.get_mag()

        results = {'numer of states': self.nstates,
                   'excitation energies': self.es_energies,
                   'excitation wavelenghts': self.es_wavelengths,
                   'oscillator strengths': self.f_osc,
                   'vel. rotatory strengths': self.r_vel,
                   'len. rotatory strengths': self.r_len,
                   'vel. transition electric dip.': self.mu_vel,
                   'len. transition electric dip.': self.mu_len,
                   'transition magnetic dip.': self.m}

        return results

    def write_stick(self):
        '''Writes Excited States data in a file for stick spectra plotting.'''

        with open('%s_stick.txt' % self.file.split('.')[0], 'w') as f:
            f.write('%5s\t%6s\t%6s\t%6s\t%10s\t%10s\n' %
                    ('State', 'E(eV)', 'l(nm)', 'f', 'R vel.', 'R len.'))

            for i in range(self.nstates):
                f.write('%5d\t%5.4f\t%5.2f\t%2.4f\t%10.4f\t%10.4f\n' %
                        (i + 1, self.es_energies[i], self.es_wavelengths[i],
                            self.f_osc[i], self.r_vel[i], self.r_len[i]))

#
# To be implemented
# I could use spectrum.py by Sandro&Lorenzo to make everything faster
# but I don't like it has too many options...I should rewrite something simpler
#  def write_spectra(self, lineshape, bandwidth):
#    '''Writes convoluted spectra in a file for plotting.'''
#    pass
#
#
# ===================
#  G09 Coupling output file
# ===================
#
# Inherits from G09 output class. This class should contain methods specific
# to extract data from G09 coupling calculations.
#
class V_output_file(output_file):
  '''A class describing G09 optimization logfiles.'''

  def __init__(self, outfile):
    output_file.__init__(self, outfile)
#
#
# ===============================
#  Generate Coupling Calculation
# ===============================
#
# This function should be called in the directory containing the directories
# of the two chromophores, which in turn contain the G09 input files of the
# eet=prop calculation
#


def gen_coup(chrom1, chrom2, opts_dict=None, thresh=15.0):
    '''Generates the coupling between chrom1 and chrom2 with options
    stored in opts_dict.'''

    f1 = input_file(os.path.join(os.getcwd(), chrom1, '%s.com' % chrom1))
    f2 = input_file(os.path.join(os.getcwd(), chrom2, '%s.com' % chrom2))

    # Check the distance between the two chromophores
    x1, y1, z1 = f1.get_com()
    x2, y2, z2 = f2.get_com()
    p1 = np.array([x1, y1, z1]) 
    p2 = np.array([x2, y2, z2])
    d = np.linalg.norm(p1-p2)

    if abs(d) > thresh:

        print("Distance %s - %s higher than %3.1f A." % (chrom1, chrom2, thresh))

    else:

        # Change the options for the generation of the new file
        opts_dict['chk'] = 'V_%s.%s' % \
            (f1.name.split('.')[0], f2.name.split('.')[0])
        
        opts_dict['title'] = 'coupling %s - %s' % \
            (f1.name.split('.')[0], f2.name.split('.')[0])
        
        opts_dict['structure'] = f1.structure + f2.structure
        
        # If the options dictionary does not refer to another coupling calculation
        # use the default options listed below for the new coupling calculation
        if 'eet=coup' not in opts_dict['job']:
            opts_dict['basis'] = f1.basis
            opts_dict['funct'] = f1.funct
            opts_dict['job'] = 'td eet=coup IOp(2/12=3) nosymm'
        
        f_coup = input_file('V_%s.%s.com' %
                            (f1.name.split('.')[0], f2.name.split('.')[0]),
                            opts_dict)
        
        # Create V_chrom1.chrom2 directory and move the .com file in it
        coup_dir = 'V_%s.%s' % (chrom1, chrom2)
        
        if not os.path.isdir(coup_dir):
            os.makedirs(coup_dir)
        
        shutil.move(f_coup.name, os.path.join(coup_dir, f_coup.name))

#
# ==============================
#  RMSD between two structures 
# ==============================
#
# The problem has been solved by W. Kabsch in 1976 (Acta Crystallographica,
# 32, 922.) The following function is described in PyMol Wiki (see
# http://www.pymolwiki.org/index.php/Kabsch). I adapted it to the structures
# obtained with this module.
#

def kabsch(struct1, struct2):
    '''Returns the RMSD calculated with Kabsch's algorithm.'''

    # Modify structures to get rid of the atomic symbol or number and convert
    # to np.array
    struct1 = np.array([ [atom[1], atom[2], atom[3]] for atom in struct1 ])    
    struct2 = np.array([ [atom[1], atom[2], atom[3]] for atom in struct2 ])    

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
    # U = np.dot(V, Wt)
 
    # rotate and translate the molecule
    # struct2 = np.dot((struct2), U)
    # struct2 = struct2 + com1

    return RMSD

#
# ============================
#  Translation of a structure 
# ============================
#

def translate(struct, dx=0.0, dy=0.0, dz=0.0):
    '''Translates a structure by dx along x, dy along y and dz along z.'''

    # Define the transformation matrix for a translation
    T = np.eye(4)
    T[-1,:3] = float(dx), float(dy), float(dz)

    # Convert the structure from the format obtained from G09 files to a format
    # for dot product with the transformation matrix
    atoms =  np.array(struct)[:,0]
    struct = np.array(struct)[:,1:].astype(np.float)

    # Add a column containing ones to the structure matrix
    struct = np.c_[struct, np.ones(len(struct))]

    # Perform the translation
    new_struct = np.dot(struct, T)

    # Return the new structure in the same format as the input one
    # I could handle the format with structured arrays, but I do not know how
    # they work, so I will do like this, for now
    final = np.c_[atoms, new_struct[:,:3]].tolist()
    final = [[atom[0], float(atom[1]), float(atom[2]), float(atom[3])] for atom
            in final]

    return final

#
# =========================
#  The Program Starts Here
# =========================
#

if __name__ == '__main__':
    pass
