#!/usr/bin/python

"""
  ++ Submit GDV ++
  - Version:     2011.11
  - Author:      Julien Bloino (idea02)
  - Description: Submission script for Gaussian
  - Edited by Lorenzo Cupellini and Sandro Jurinovich

"""

import sys, os, re
from subprocess import Popen, PIPE

# ==================================================
#  MACHINE-SPECIFIC VARIABLES
# ==================================================

SCRATCHDIR  = '/local/'
QUEUES      = {'short_old' : (12,  '2500MW'),
               'long_old'  : (12,  '2500MW'),
               'short_new' : (16,  '7000MW'),
               'long_new'  : (16,  '7000MW'),
               'raptor'    : (12,  '128GB'),
               'gpu'       : (16,  '7000MW')}

# ==================================================
#  GLOBAL VARIABLES
# ==================================================

jobPID      = str(os.getpid())
PROGNAME    = os.path.basename(sys.argv[0])
USER        = os.getenv('USER')
SCRATCHDIR  = os.path.join(SCRATCHDIR,USER,'gaurun'+jobPID)
STARTDIR    = os.getcwd()
GAUSSIANDIR = '/cm/shared/apps/gaussian'

# ==================================================
#  FUNCTIONS
# ==================================================

def GetOpts():
  """
  Analyzes the command line and sets adequately the runtime variables
  It also retrieves the filenames from the arguments
  """
  subOpts = {'queue': None,  # Queue name for job submission
             'job'  : None,  # Job name
             'quiet': True, # Silent mode
             'add'  : None,  # Add qsub command lines
             'expert': 0}    # Expert mode for running Gaussian
  gauOpts = {'root'  : None, # Root dir. where the Gaussian version is stored
             'exe'   : None, # Gaussian executable (the alias for Link0)
             'out'   : None, # Filename for the output of Gaussian
             'chk'   : None, # Filename for the Gaussian checkpoint
             'keep'  : {'chk' : False,  # Keep chk filename given in input file
                        'mem' : False,  # Keep mem request given in input file
                        'proc': False}, # Keep num. of proc. from input file
             'force' : {'chk': False,   # control that chk file exists
                        'FC' : False},  # control that files exist for FC calc.
             'ignore': {'chk': False},  # ignore checkpoint file in all treatment
             'wrkdir': None} # Working directory (for dev-related works)

  gauRev = {'g03c01'     : os.path.join(GAUSSIANDIR,'g03C01','g03'),
            'g09d01'     : os.path.join(GAUSSIANDIR,'g09D01','g09'),
            'gdvh11'     : os.path.join(GAUSSIANDIR,'gdvH11','gdv'),
            'gdvh11eet'  : os.path.join(GAUSSIANDIR,'gdvH11eet','gdv'),
            'gdvh23'     : os.path.join(GAUSSIANDIR,'gdvH23','gdv'),
            'gdvh23dev'  : os.path.join(GAUSSIANDIR,'gdvH23dev','gdv'),
            'gdvh28'     : os.path.join(GAUSSIANDIR,'gdvH28','gdv'),
            'gdvh33p'    : os.path.join(GAUSSIANDIR,'gdvH33p','gdv'),
            'gdvh33p_old': os.path.join(GAUSSIANDIR,'gdvH33p_acml_nehalem','gdv'),
            'gdvh33p_new': os.path.join(GAUSSIANDIR,'gdvH33p_acml_sandy','gdv'),
            'gdvh33p_dev': os.path.join(GAUSSIANDIR,'gdvH33p_dev','gdv'),
            'gdvh36eet': os.path.join(GAUSSIANDIR,'gdvH36eet','gdv'),
            'gdvh36m': os.path.join(GAUSSIANDIR,'gdvH36-molecolab','gdv'),
            'gdvh36': os.path.join(GAUSSIANDIR,'gdvH36','gdv')}

  # Check if help has been requested
  if set(['-h', '--help']) & set(sys.argv):
    print """Usage: %s [options] <Gaussian_input_file>
 Options are:
  ## QSUB-RELATED OPTIONS ##
  -q, --queue:
    Sets the queue type. The possible choices are:
    + enlight    ( 3 days,    12 cores,  24GB RAM)
    + long       ( unlimited, 12 cores,  24GB RAM)
    + big        ( unlimited, 12 cores, 144GB RAM)

  -j, --job:
     Sets the job name.

  -s, --silent:
     Do not save standard output and error in files
     WARNING: The consequence will be a loss of these outputs

  !! Expert use !!
  -0:
     Runs input without modifying or analyzing it
     This option deactivates all controls and simply run Gaussian with the given input.
     Parameters for Link0 should be correctly given by the user.

  ## GAUSSIAN-RELATED OPTIONS ##
  -g, --gauroot:
     Sets the path to the Gaussian executables. It can be given as an absolute
     path or with the following keywords:
     + g03c01    : Gaussian 03 Rev. C.01 
     + g09d01    : Gaussian 09 Rev. D.01 (default)
     + gdvh11    : Gaussian DV Rev. H.11
     + gdvh11eet : Gaussian DV Rev. H.11 + EET
     + gdvh23    : Gaussian DV Rev. H.23
     + gdvh23der : Gaussian DV Rev. H.23 + DER
     + gdvh28    : Gaussian DV Rev. H.28
     + g09       : alias for g09d01
     + gdv       : alias for gdvh28

  -o, --out:
     Sets the output filename.

  -c, --chk:
     Sets the checkpoint filename.

  -k, --keep:
     Keeps user-given parameters to control the Gaussian job in input file.
     The possible options are:
     + c, chk: Keeps checkpoint filename given in input
     + m, mem: Keeps memory requirements given in input
     + p, proc: Keeps number of proc. required in input
     + a, all: Keeps all data liste above
     NOTE: Several options can be given separated by commas.

  -f, --force:
     Force controls with respect to the data in Gaussian input:
     + c, chk: if checkpoint file does not exist, aborts
     + FC: if input files for FC calc do not exist, aborts
     NOTE: This is useful to check that the calculations will not abort
           before the actual calculations

  -i, --ignore:
     Ignore the following options in the command list:
     + c, chk: ignore the checkpoint file (omit it in the input, do not copy it)

  -w, --wrkdir:
     Appends a working directory to the Gaussian path to look for executables.
     Several working directories can be given by using multiple times the -w
       option. In this case, the last given working directory is the first one
       using during the path search.
     NOTE: Only the path to the working root directory is needed. The script
           automatically selects the correct directories from there.
     WARNING: The script expects a working tree structure as intended in the
              Gaussian developer manual.

  NOTES:
   - Suboptions must be given separated by commas without spaces
   - Short (-s) and long (--long) options follow the typical UNIX behavior:
     + -s value OR -svalue
     + --long=value

  EXAMPLES OF USAGE:
   - Direct submission, using default parameters
     %s input.com
   - Development-related tests on a short queue
     %s -q cray_short -w /home/idea02/working -g gdvh09 input.com
   - Keep memory and processors-related information from input file
     %s -km,p input.com
   - Run FC-related jobs controlling the source files are all available
     %s -fFC -kc input.com
     Note: For typical FC calculations, -fFC == -fFC -kc
   - Explicit path to the Gaussian executable
     %s -g /share/gaussian/g09.b01/g09 input.com
""" % tuple([PROGNAME for i in range(6)])
    sys.exit()

  i = 1
  lenArgs = len(sys.argv)
  queue   = ''
  job     = ''
  gauroot = ''
  outfile = ''
  chkfile = ''
  keep    = None
  force   = None
  ignore  = None
  wrkdir  = []
  files   = []

  while i < lenArgs:
    # OPTIONS STARTING WITH '--'
    if sys.argv[i].startswith('--'):
      arg = sys.argv[i][2:]
      if arg.find('=') <= 0:
        if arg == 'silent': subOpts['quiet'] = True
        else:
          print 'ERROR: Incorrect form for argument: "%s"' % sys.argv[i]
          sys.exit()
      else:
        option, value = sys.argv[i][2:].split('=')
        if not value:
          if   option == 'keep' : keep = True
          elif option == 'force': force = True
          else:
            print 'ERROR: A value is expected: --%s=value' % option
            print 'See help for the possible values'
            sys.exit()
        else:
          if   option == 'queue'  : queue   = value
          elif option == 'job'    : job     = value
          elif option == 'gauroot': gauroot = value
          elif option == 'out'    : outfile = value
          elif option == 'chk'    : chkfile = value
          elif option == 'keep'   : keep    = value
          elif option == 'force'  : force   = value
          elif option == 'ignore' : ignore  = value
          elif option == 'wrkdir' : wrkdir.append(value)
          else:
            print 'Option "'+option+'" is not recognized.'
            sys.exit()

    # OPTIONS STARTING WITH '-'
    elif sys.argv[i].startswith('-'):
      option = sys.argv[i][1]
      value = ''
      # Options without values are treated first
      if option == 's':
        subOpts['quiet'] = True
      elif option == '0':
        subOpts['expert'] = 1
      # Remaining options with possible values
      else:
        if len(sys.argv[i]) > 2: # Value is included in arg
          value = sys.argv[i][2:]
        else:
          if i < lenArgs-1 and not sys.argv[i+1].startswith('-'):
            # There is a possible value only if Option is not the last argument
            # and the following argument does not start with a dash
            # NOTE: for now, we discard the possibility of negative number as val
            value = sys.argv[i+1]
            if option in ['f','k'] and os.path.exists(value):
              value = None
            else:
              i += 1
  
        if not value:
          if option == 'f':   force = True
          elif option == 'k': keep  = True
          else:
            print 'Value is missing for option "'+option+'"'
            sys.exit()
        else:
          if   option == 'c': chkfile = value
          elif option == 'f': force   = value
          elif option == 'g': gauroot = value
          elif option == 'i': ignore  = value
          elif option == 'j': job     = value
          elif option == 'k': keep    = value
          elif option == 'o': outfile = value
          elif option == 'q': queue   = value
          elif option == 'w': wrkdir.append(value)
          else:
            print 'Option "'+option+'" is not recognized.'
            sys.exit()

    # OTHER ARGUMENTS
    else:
      files.append(sys.argv[i])
    # Increment to following argument
    i += 1

  ### ANALYSIS OF THE OPTIONS ###
  ###############################
  if queue:
    if queue in QUEUES.keys(): subOpts['queue'] = queue
    else:
      print 'A queue identifier is required. See help for details.'
      sys.exit()
  else:
    subOpts['queue'] = 'short_new'

  if gauroot:
    if gauroot in gauRev.keys(): gauOpts['root'] = gauRev[gauroot]
    elif gauroot == 'gdv': gauOpts['root'] = gauRev['gdvh28']
    elif gauroot == 'g09': gauOpts['root'] = gauRev['g09d01']
    else:
      if os.path.exists(gauroot): gauOpts['root'] = gauroot
      else:
        print 'Path "'+gauroot+'" does not exist'
        sys.exit()
  else:
    gauOpts['root'] = gauRev['g09d01']

  basedir, gauOpts['exe'] = os.path.split(gauOpts['root'])
  # Remove possible duplicate if executable has been given instead of root
  if os.path.split(basedir)[-1] == gauOpts['exe']:
    gauOpts['root'] = os.path.split(gauOpts['root'])[0]
  if gauOpts['exe'] not in ['g03', 'gdv', 'g09']:
    print 'Gaussian version is not recognized: must be g03, g09 or gdv'
    sys.exit()

  if force:
    if force == True:
      gauOpts['force']['chk'] = True
      gauOpts['force']['FC']  = True
    else:
      values = force.split(',')
      for value in values:
        if value in ['c', 'chk']: gauOpts['force']['chk'] = True
        elif value == 'FC': gauOpts['force']['FC'] = True
        else:
          print 'Unrecognized choice for force:', value
          sys.exit()
      
  if keep:
    if keep == True:
      gauOpts['keep']['chk']  = True
      gauOpts['keep']['mem']  = True
      gauOpts['keep']['proc'] = True
    else:
      values = keep.split(',')
      for value in values:
        if value in ['c', 'chk']:    gauOpts['keep']['chk'] = True
        elif value in ['m', 'mem']:  gauOpts['keep']['mem'] = True
        elif value in ['p', 'proc']: gauOpts['keep']['proc'] = True
        elif value in ['a', 'all']:
          gauOpts['keep']['chk']  = True
          gauOpts['keep']['mem']  = True
          gauOpts['keep']['proc'] = True
        else:
          print 'Unrecognized choice for keep:', value
          sys.exit()

  if ignore:
    if ignore == True:
      gauOpts['ignore']['chk']  = True
    else:
      values = ignore.split(',')
      for value in values:
        if value in ['c', 'chk']: gauOpts['ignore']['chk'] = True

  if wrkdir:
    gauOpts['wrkdir'] = []
    for value in wrkdir:
      if os.path.exists(value): gauOpts['wrkdir'].append(value)
      else:
        print 'Working directory "'+value+'" does not exist'
        sys.exit()

  # We need to identify the input file for the remaining options
  if len(files) > 1:
    print 'Only one input file must be given'
    sys.exit()
  
  gauInp = files[0]
  # Check that Gaussian Input file exists
  if not os.path.exists(gauInp):
    print 'Error: Cannot find Gaussian input file "'+gauInp+'"'
    sys.exit()
  filename = os.path.basename(gauInp)
  basename = os.path.splitext(filename)[0]

  if subOpts['expert'] != 1:
    if chkfile:
      gauOpts['chk'] = chkfile
    elif not gauOpts['keep']['chk']:
      gauOpts['chk'] = basename + '.chk'

  if outfile: gauOpts['out'] = outfile
  else: gauOpts['out'] = basename + '.log'

  if job: subOpts['job'] = job
  else: subOpts['job'] = basename
       
  return gauInp, subOpts, gauOpts

# ----------------------------------------

def chkGauInp(refInp,newInp,NProcs,Mem,Chk,DoCopy,Omit):
  io1 = open(refInp, 'r')
  io2 = open(newInp, 'w')

  gauDat = {'NProcs': NProcs, 'Mem': Mem, 'Chk': Chk}
  subDat = []

  L1ToRead  = True
  FCInputs  = ['chk', 'inp', 'rwf']
  FCFiles   = []
  fmtFrq    = r'\bfreq\w*=?(?P<delim>\()?\S*'+'%s'+r'\S*(?(delim)\))\b'
  fmtGeom   = r'\bgeom\w*=?(?P<delim>\()?\S*'+'%s'+r'\S*(?(delim)\))\b'
  strFC     = r'\b(fc|fcht|ht)\b'
  keyAllChk = re.compile(fmtGeom % r'\ballcheck\b')
  keyGeom   = re.compile(fmtGeom % r'\bcheck\b')
  keyFC     = re.compile(strFC)
  keyFrqFC  = re.compile(fmtFrq % strFC)
  keyFCOpt  = re.compile(fmtFrq % r'\breadfcht\b')
  keyAnhOpt = re.compile(fmtFrq % r'\breadanh')
  keyAnharm = re.compile(fmtFrq % r'\banharm(|onic)\b')
  keyFCInp  = re.compile(r'\b(chk1|chk2|out1|out2)\b')
  FCArgs    = r'adiabaticfc|adiabaticshift|verticalfc|lcm|verticalgradient|'   \
            + r'fc|fcht|ht|op|one-photon|ecd|abs|absorption|emi|emission|'     \
            + r'calc1|chk1|out1|calc2|chk2|out2|nstate|inpfrq1|inpfrq2|'       \
            + r'noreaddip|inpedip|inpmdip|inpdipa|inpdipb|sclvec|inpdener|'    \
            + r'refstate|jdusch|jident|jreduced|jdiag|blockthresh|blocktol|'   \
            + r'maxc1|maxc2|maxovr|maxcmb|noreli00|maxint|gaussian|'           \
            + r'lorentzian|stickspectrum|specmin|specmax|specres|spechwhm|'    \
            + r'allspectra|prtspectra|prtmat|prtint|dotemp|minpop|'            \
            + r'temperature|maxstatesi|transposedht|nointan|maxband1|maxosc1|' \
            + r'maxband2|maxosc2|maxbands|deltasp|clearanhvfc|forcefccalc|'    \
            + r'forceprtspectrum|reorient|noreorient|rotation|rotniter'
  keyFCArgs = re.compile(r'\b('+FCArgs+r')\b')
  keyGeomView = re.compile(r'\bgeomview\b')
  DoFCHT    = False
  ReadFCHT  = False
  DoAnharm  = False
  ReadAnh   = False
  block     = ''
  nblock    = 0
  
  # Read Input File
  if Mem:    io2.write('%%Mem=%s\n' % Mem)
  if NProcs: io2.write('%%NProcShared=%s\n' % NProcs)
  if Chk:    io2.write('%%Chk=%s\n' % Chk)
  for Line in io1:
    Line2 = Line.strip().lower()
    # ------------------------------
    # Reset on new Gaussian job
    #
    if Line2 == '--link1--':
      L1ToRead  = True
      DoFCHT    = False
      DoAnharm  = False
      ReadFCHT  = False
      ReadAnh   = False
      if Mem:    Line += '%%Mem=%s\n' % Mem
      if NProcs: Line += '%%NProcShared=%s\n' % NProcs
      if Chk and not Omit['chk']:
        Line += '%%Chk=%s\n' % Chk
    # ------------------------------
    # Link 0 parameters
    #
    elif Line2.startswith('%'):
      if Line2.startswith('%chk'):
        if not Chk: gauDat['Chk'] = Line.split('=')[1].strip()
        elif not Omit['chk']: Line = ''
      elif Line2.startswith('%mem'):
        if Mem: Line = ''
        else: gauDat['Mem'] = Line.split('=')[1].strip()
      elif Line2.startswith('%nproc'):
        if NProcs: Line = ''
        else: gauDat['NProcs'] = int(Line.split('=')[1].strip())
    # ------------------------------
    # Route parameters (Link 1)
    #
    elif Line2.startswith('#') and L1ToRead:
      L1ToRead = False
      while Line.strip() != '':
        io2.write(Line)
        Line = io1.next()
        Line2 += Line.strip().lower()
      if re.search(keyFCOpt, Line2): ReadFCHT = True
      if re.search(keyFrqFC, Line2): DoFCHT = True
      if re.search(keyAllChk, Line2):
        nblock = 0
      elif re.search(keyGeom, Line2):
        nblock = 2
      if ReadFCHT and not DoFCHT:
        print 'ERROR: ReadFCHT is not sufficient to run FCHT calculations'
        io1.close()
        io2.close()
        sys.exit()
      if re.search(keyAnharm, Line2): DoAnharm = True
      if re.search(keyAnhOpt, Line2): ReadAnh = True
      if ReadAnh and not DoAnharm:
        print 'ERROR: ReadAnharm is not sufficient to run anharmonic calc.'
        sys.exit()
      if DoAnharm and DoFCHT:
        print 'SORRY, Simultaneous treatment of FC and Anharm is not yet ' \
              + 'available. Please contact the author if needed.'
        sys.exit()
    # ------------------------------
    # Screen input after route
    #
    #Julien: we could add treatment for @ if we want to check file exists
    elif Line2 and not L1ToRead:
      block = ''
      while Line.strip() != '':
        io2.write(Line)
        block += Line
        Line = io1.next()
        Line2 += Line.strip().lower()
      if re.search(keyGeomView, Line2):
        subDat.append(['cpfrom', 'points.off'])
      if ReadFCHT and nblock == 0:
        if re.search(keyFCArgs, Line2):
          io2.write(Line)
          Line = io1.next()
          while Line.strip() != '':
            FCFiles.append(Line.strip())
            io2.write(Line)
            Line = io1.next()
      elif nblock > 0:
        nblock -= 1
    io2.write(Line)

  # We assume that the last block contains the files to copy
  if DoFCHT and not ReadFCHT:
    for line in block.split('\n'):
      if line.strip(): FCFiles.append(line.strip())

  io1.close()
  io2.close()

  # Copy files for CHK
  if not Omit['chk']:
    if not os.path.exists(gauDat['Chk']):
      if DoCopy['chk'] or (DoFCHT and len(FCFiles)==1 and DoCopy['FC']):
        print 'Error: Checkpoint file "'+ gauDat['Chk']+'" does not exist'
    else:
      subDat.append(['cpto', gauDat['Chk']])
 
  #lore
  # Copy EET files
  EETFiles = ['CisPATr','CisPATr2','CisPA','CisPA2',
              'CisPBTr','CisPBTr2','CisPB','CisPB2',
             'Dat','Dat2','SelCoup.dat']

  for EETFile in EETFiles:
    if os.path.exists(EETFile):
      subDat.append(['cpto', EETFile])

  for EETFile in EETFiles:
    if os.path.exists(EETFile):
      subDat.append(['cpfrom', EETFile])

  # FragData.dat
  #subDat.append(['cpfrom','CisAttA'])
  #subDat.append(['cpfrom','CisDetA'])
  #subDat.append(['cpfrom','FragData.dat'])
  #subDat.append(['cpfrom','CisPATr'])
  #subDat.append(['cpfrom','CisPA'])
  #subDat.append(['cpfrom','CisPBTr'])
  #subDat.append(['cpfrom','CisPB'])
  #subDat.append(['cpfrom','Dat'])

  # Franck-Condon calculations: first check everything is ok with FCFiles
  if DoFCHT and len(FCFiles) == 0:
    print 'ERROR: I Could not find files for Franck-Condon calculations'
    print '       I will abort now.'
    sys.exit()
  # copy files for FCHT
  for FCFile in FCFiles:
    if not os.path.exists(FCFile):
      if DoCopy['FC']:
        print 'ERROR: The following file is missing to carry out FC calc.:',
        print FCFile
        sys.exit()
    else:
      subDat.append(['cpto', FCFile])  

  return gauDat, subDat

# ----------------------------------------

def setGauEnv(GauDir):
  l = [os.path.join(GauDir,f) for f in ['bsd','local','extras','']]
  GAUSS_EXEDIR = os.pathsep.join(l)
  TORQUE_DIR   = os.pathsep.join(['/cm/shared/apps/torque/4.2.2/bin','/cm/shared/apps/torque/4.2.2/sbin'])
  os.environ['GAUSS_EXEDIR'] = GAUSS_EXEDIR
  os.environ['GAU_ARCHDIR'] = os.path.join(GauDir,'arch')
  os.environ['PATH'] = os.pathsep.join([GAUSS_EXEDIR,TORQUE_DIR]) #os.pathsep.join([GAUSS_EXEDIR,os.environ['PATH']])
# os.environ['LD_LIBRARY_PATH'] = GAUSS_EXEDIR #\
#    os.pathsep.join([GAUSS_EXEDIR,os.environ['LD_LIBRARY_PATH']])

# ==================================================
#  MAIN PROGRAM
# ==================================================
if __name__ == '__main__':
  if len(sys.argv) < 2:
    print 'Error: Input file is missing.'
    print 'Type "'+PROGNAME+' -h" for more information'
    sys.exit()

  gauInp, subOpts, gauOpts = GetOpts()
  newDir, refInp = os.path.split(gauInp)
  WORKDIR = os.path.abspath(newDir)

  # A new, temporary input is created
  newInp = '%s_%s.com' % (os.path.splitext(refInp)[0], jobPID)
  basename = os.path.splitext(refInp)[0]

  # Define NProcs and Mem
  NProcs = None
  Mem    = None
  if not gauOpts['keep']['proc']:
    NProcs = QUEUES[subOpts['queue']][0]
  if not gauOpts['keep']['mem']:
    Mem = QUEUES[subOpts['queue']][1]
    
  # The script works in the directory where the input file is stored
  os.chdir(WORKDIR)
  if subOpts['expert'] == 0:
    gauDat, subDat = chkGauInp(refInp,newInp,NProcs,Mem,gauOpts['chk'],
      gauOpts['force'],gauOpts['ignore'])
    if not gauOpts['chk']: gauOpts['chk'] = gauDat['Chk']

    if gauDat['NProcs'] > QUEUES[subOpts['queue']][0]:
      print 'Too many processors required for the chosen queue'
      sys.exit()
  elif subOpts['expert'] == 1:
    gauDat['NProcs'] = QUEUES[subOpts['queue']][0]

  # All data have been gathered for calculation. Set environment variables
  setGauEnv(gauOpts['root'])
  
  # ==================================================
  #  PBS-RELATIVE VARIABLES
  # ==================================================
  # NOTA: For the sake of simplicity, the script SHELL is forced to be TCSH
  
  PBSHEADER = """
echo "----------------------------------------"
echo "PBS queue:     "$PBS_O_QUEUE
echo "PBS host:      "$PBS_O_HOST
echo "PBS node:      "$HOSTNAME
echo "PBS workdir:   %s"
echo "PBS jobid:     "$PBS_JOBID
echo "PBS jobname:   "$PBS_JOBNAME
echo "PBS inputfile: %s"
echo "PATH:          "$PATH
echo "LD_LIB_PATH    "$LD_LIBRARY_PATH
echo "gaucmd         "
echo "GAUSS_EXEDIR:  "$GAUSS_EXEDIR
echo "----------------------------------------"
  """ % (SCRATCHDIR,gauInp)
  
  PBSCMDS = 'mkdir -p %s\n' % SCRATCHDIR
  PBSCMDS += 'cd %s\n' % SCRATCHDIR
  # Move temporary input file
  PBSCMDS += 'mv %s ./\n' % os.path.join(WORKDIR,newInp)
  # Copy files given in Input file if available
  if len(subDat) > 0: 
    for data in subDat:
      if (data[0] == 'cpto'):
        PBSCMDS += '(cp %s ./) >& /dev/null\n' % os.path.join(WORKDIR,data[1])

  gauArgs=''
  if gauOpts['wrkdir']:
    if len(gauOpts['wrkdir']) > 0:
      gauArgs += ' -exedir="'
      for value in gauOpts['wrkdir']:
        if os.path.exists(value): gauArgs += '%s/l1:%s/exe-dir:' % \
          (value, value)
      gauArgs += '$GAUSS_EXEDIR"'

  PBSCMDS += '%s %s < %s > %s\n' % (gauOpts['exe'],gauArgs.strip(),newInp,
    os.path.join(WORKDIR,gauOpts['out']))
      
  if not gauOpts['ignore']['chk']:
    PBSCMDS += '(cp %s %s) >& /dev/null\n' % (gauOpts['chk'], WORKDIR)
  if len(subDat) > 0: 
    for data in subDat:
      if (data[0] == 'cpfrom'):
        if (data[1] == 'FragData.dat'):
          PBSCMDS += '(cp %s %s) >& /dev/null\n' % (data[1], os.path.join(WORKDIR,basename+'.dat'))
        else:
          PBSCMDS += '(cp %s %s) >& /dev/null\n' % (data[1], WORKDIR)
  
  # Cleaning
  PBSCMDS += 'cd .. \nrm -rf %s\n' % SCRATCHDIR
  
  # ==================================================
  #  SUBMISSION JOB
  # ==================================================
  subCmd = 'qsub' # QSub is used to send the job
  subCmd += ' -V' # exports user environment variables
  subCmd += ' -r n' # Flags the current job as non-rerunnable if it fails due
                    # to an error of the working node (safer than force rerun)
  subCmd += ' -N \'%s\'' % subOpts['job'].replace(' ','_') # Sets jobname
  if subOpts['queue'] == 'short_old':
    subCmd += ' -l nodes=1:ppn=12:oldres'
  elif subOpts['queue'] == 'long_old':
    subCmd += ' -l nodes=1:ppn=12:oldres'
  elif subOpts['queue'] == 'short_new':
    subCmd += ' -l nodes=1:ppn=16:newres'
  elif subOpts['queue'] == 'long_new':
    subCmd += ' -l nodes=1:ppn=16:newres'
  elif subOpts['queue'] == 'raptor':
    subCmd += ' -l nodes=1:ppn=12:raptorres'
  elif subOpts['queue'] == 'gpu':
    subCmd += ' -l nodes=1:ppn=16:gpures'
    # Nodes-related option: fixes the number of nodes to use (nodes=), and the
    # number of processors per node (ppn=)
  subCmd += ' -q %s' % subOpts['queue'] # Queue name
  if (subOpts['quiet']): # Quiet mode, all outputs are redirected to /dev/null
    subCmd += ' -o /dev/null -e /dev/null'
#  subCmd += ' -' # Ends QSub options

  PBSHEADER += 'echo "\n   === LIST OF COMMANDS ===\n"\necho "' \
    + PBSCMDS.replace('\n','"\necho "') + '"\n'
# print subCmd
# print PBSHEADER
# print PBSCMDS
  subSend = Popen(args=subCmd, shell=True, stdin=PIPE, stdout=PIPE)
  subSend.stdin.write(PBSHEADER)
  subSend.stdin.write(PBSCMDS)
  subSend.stdin.close()
  print "Comando: "+subCmd
  print "QSub submission job : '"+subSend.stdout.readline().strip()+"'"

