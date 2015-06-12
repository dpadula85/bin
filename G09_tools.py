#!/usr/bin/env python

# Import standard Python modules
import sys
import os
import re
import shutil
import argparse as arg

# Import personal modules
import G09_files


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Processes G09 files.')

    # Positional arguments
    # The same job will be performed on all these files
    # The argument is mandatory, but some following options (e.g. --coup)
    # do not require it
    parser.add_argument('f_inp', nargs='+', help='''Names of the G09 logfiles.
    Some options do not require it to be specified.''')

    # Optional arguments
    # Extraction of the Optimized Geometry and creation of a new G09 input file
    parser.add_argument('--opt', help='''Extracts the optimized geometry and
    creates a new G09 input file with it.''', action='store_true')

    # Extraction of Excited States properties and creation of a stick spectrum
    parser.add_argument('--stick', help='''Extracts Excited States
    properties and writes it to a stick file.''', action='store_true')

    # Preparation of input files for Coupling calculations for EXAT
    parser.add_argument('--coup', help='''Prepares Coupling files and directories
    for coupling calculations with G09.
    This option does not require a file to be specified.''',
                        choices=['all', 'calclist'])

    # Preparation of EXAT calculations including EXAT calculations for each
    # couple of chromophores constituting the system, or not.
    # bb is for a backbone calculation, including only n->pi* and pi->pi*
    # transitions described by states 1 and 4 of NMA.
    parser.add_argument('--exat', help='''Prepares files needed for EXAT.
    This option does not require a file to be specified.''',
                        choices=['couples', 'system', 'bb'])

    # Section of options for G09 - ***DEFAULT OPTIONS ARE SET HERE***
    # This will be called in the main program
    parser.add_argument('-n', '--nproc', default=16, type=int,
                        help='Number of processors.')

    parser.add_argument('-m', '--mem', default=24, type=int,
                        help='Memory in GB.')

    parser.add_argument('-chk', '--chk',
                        help='Name of the checkpoint file without extension.')

    parser.add_argument('-j', '--job', nargs='*', default=['td=(nstates=32)'],
                        required=False, help='''Job type. This does not work combined
    with --coup option.''')

    parser.add_argument('-f', '--funct', default='cam-b3lyp',
                        help='Functional.')

    parser.add_argument('-b', '--basis', default='6-31+G(d)',
                        help='Basis set.')

    parser.add_argument('-c', '--charge', default=0, type=int,
                        help='Charge.')

    parser.add_argument('-mt', '--mult', default=1, type=int,
                        help='Multiplicity.')

    # Print help if no arguments are received
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    args.job = ' '.join(args.job)

    G09_opts = {'name': '',
                'nproc': args.nproc,
                'mem': args.mem,
                'mem_unit': 'GB',
                'chk': args.chk,
                'funct': args.funct,
                'basis': args.basis,
                'job': args.job,
                'title': 'Title card',
                'charge': args.charge,
                'mult': args.mult,
                'structure': [],
                'add_opts': '\n'}

    return args, G09_opts

#
# =========================
#  The Program Starts Here
# =========================
#

if __name__ == '__main__':

    # Parse the options
    args, G09_opts = options()

    # Analyze each option
    if args.opt:

        for f_inp in args.f_inp:

            f = G09_files.opt_output_file(f_inp)
            G09_opts['structure'] = f.structure
            G09_opts['chk'] = '%s_new.com' % f.name.split('.')[0]
            G09_files.input_file('%s_new.com' % f.name.split('.')[0], G09_opts)

    if args.stick:

        for f_inp in args.f_inp:

            f = G09_files.td_output_file(f_inp)
            f.write_stick()

    if args.coup:

        # Directories of chromophores, expressed as all directories whose name
        # contains only digits
        dirs = sorted(filter(
            lambda x: re.match('\d+', x), os.listdir(os.getcwd())))

        if args.coup == 'all':

            # Write calclist.in file for properties and coupling calculations
            calclist = open('calclist.in', 'w')
            calclist.write('1 1 0 0\n')

            G09_opts['funct'] = 'b3lyp'
            G09_opts['job'] = '''
            integral(ultrafine) td eet=coup IOp(2/12=3) nosymm'''

            for chrom1 in dirs[:]:
                for chrom2 in filter(lambda x: x != chrom1, dirs):

                    G09_files.gen_coup(chrom1, chrom2, G09_opts)

                # Remove chrom1 to avoid equivalent couple chrom2.chrom1
                calclist.write('RES: %s\n' % chrom1)
                dirs.remove(chrom1)

            calclist.close()

        if args.coup == 'calclist':

            calclist = open('calclist.in', 'w')
            calclist.write('1 1 0 0\n')

            for chrom in dirs:
                calclist.write('RES: %s\n' % chrom)

            calclist.close()

    # exat option has three possibilities
    if args.exat:

        # Directories of chromophores, expressed as all directories whose name
        # contains only digits
        dirs = sorted(filter(
            lambda x: re.match('\d+', x), os.listdir(os.getcwd())))

        if args.exat == "system":

            chromlist = open('chromlist.in', 'w')
            for chrom in dirs:
                chromlist.write('%s\n' % chrom)

            chromlist.close()

        if args.exat == "bb":

            chromlist = open('chromlist.in', 'w')
            for chrom in dirs:
                chromlist.write('%s 1 4\n' % chrom)

            chromlist.close()

        if args.exat == "couples":

            chromlist = open('chromlist.in', 'w')

            for chrom1 in dirs[:]:
                for chrom2 in filter(lambda x: x != chrom1, dirs):

                    # Create a directory called exat_chrom1_chrom2 and copy
                    # chrom1, chrom2 and V_chrom1.chrom2 directories
                    exat_dir = 'exat_%s.%s' % (chrom1, chrom2)
                    coup_dir = 'V_%s.%s' % (chrom1, chrom2)
                    os.makedirs(os.path.join(exat_dir, chrom1))
                    os.makedirs(os.path.join(exat_dir, chrom2))
                    os.makedirs(os.path.join(exat_dir, coup_dir))

                    shutil.copyfile(os.path.join(chrom1, '%s.log' % chrom1),
                                    os.path.join(exat_dir, chrom1, '%s.log' %
                                                 chrom1))

                    shutil.copyfile(os.path.join(chrom2, '%s.log' % chrom2),
                                    os.path.join(exat_dir, chrom2, '%s.log' %
                                                 chrom2))

                    shutil.copyfile(os.path.join(coup_dir, '%s.log' % coup_dir),
                                    os.path.join(exat_dir, coup_dir, '%s.log' %
                                                 coup_dir))

                    # Write chromlist.in in the exat_dir for the partial job
                    chromlist_part = open(os.path.join(
                        exat_dir, 'chromlist.in'), 'w')

                    chromlist_part.write('%s\n' % chrom1)
                    chromlist_part.write('%s\n' % chrom2)
                    chromlist_part.close()

                # Remove chrom1 to avoid equivalent couple chrom2.chrom1
                chromlist.write('%s\n' % chrom1)
                dirs.remove(chrom1)

            chromlist.close()
