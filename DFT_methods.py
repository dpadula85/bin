#!/usr/bin/env python

#DFT Functionals in G09
#From http://www.gaussian.com/g_tech/g_ur/k_dft.htm


functionals = {'exchange'     : ['S', 'XA', 'B', 'PW91', 'mPW', 'G96',
                               'PBE', 'O', 'TPSS', 'BRx', 'PKZB', 'wPBEh',
                               'PBEh'],
               'correlation'  : ['VWN', 'VWN5', 'LYP', 'PL', 'P86', 'PW91',
                               'B95', 'PBE', 'TPSS', 'KCIS', 'BRC', 'PKZB',
                               'VP86', 'V5LYP'],
               'only_xc'      : ['HFS', 'XAlpha', 'HFB'],
               'pure'         : ['VSXC', 'HCTH', 'HCTH93', 'HCTH147', 'HCTH407',
                                  'tHCTH', 'M06L', 'B97D', 'B97D3', 'SOGGA11', 'M11L',
                                  'N12', 'MN12L'],
               'hybrid'       : ['B3LYP', 'B3P86', 'B3PW91', 'B1B95', 'mPW1PW91',
                                 'mPW1LYP', 'mPW1PBE', 'mPW3PBE', 'B98', 'B971', 'B972',
                                 'PBE1PBE', 'B1LYP', 'O3LYP', 'BHandH', 'BHandHLYP', 'BMK',
                                 'M06', 'M06HF', 'M062X', 'tHCTHhyb', 'APFD', 'APF', 'SOGGA11X',
                                 'PBEh1PBE', 'TPSSh', 'X3LYP'],
               'long_range'   : ['HSEH1PBE', 'OHSE2PBE', 'OHSE1PBE', 'wB97XD', 'wB97', 'wB97X',
                                 'LC-wPBE', 'CAM-B3LYP', 'HISSbPBE', 'M11', 'N12SX', 'MN12SX']}


functionals['xc_corr'] = [ xc + corr for xc in functionals['exchange'] for corr in functionals['correlation']]


#Basis sets in G09
#From http://www.gaussian.com/g_tech/g_ur/m_basis_sets.htm


basis = {'pople'    : ['STO-3G', '3-21G', '6-21G', '4-31G', '6-31G', '6-311G'],
         'dunnings' : ['cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ','cc-pV5Z', 'cc-pV6Z'],
         'alrichs'  : ['SV', 'SVP', 'TZV', 'TZVP', 'QZVP'],
         'other'    : ['D95', 'D95V', 'SHC', 'CEP-4G', 'CEP-31G', 'CEP-121G', 'LanL2MB',
                       'LanL2DZ', 'SDD', 'SDDAll', 'Def2','MidiX', 'EPR-II', 'EPR-III',
                       'UGBS', 'MTSmall', 'DGDZVP', 'DGDZVP2','DGTZVP', 'CBSB7']}


#Dictionary of job types
#From http://www.gaussian.com/g_tech/g_ur/m_jobtypes.htm


job_types = {'opt'              : 'opt',
             'opt=modredundant' : 'scan',
             'scan'             : 'scan',
             'freq'             : 'freq',
             'vcd'              : 'vcd',
             'roa'              : 'roa',
             'td'               : 'td',
             'eet=prop'         : 'prop',
             'eet=coup'         : 'coup',
             'force'            : 'gradients',
             'optrot'           : 'optical rotation'}


#=========================
# The Program Starts Here
#=========================

if __name__ == '__main__':
  
  print
  print '''DFT Functionals implemented in G09.\nSource: http://www.gaussian.com/g_tech/g_ur/k_dft.htm'''
  print
  for key, value in functionals.iteritems():
    print '%s : \n' % key
    for funct in value:
      print funct,
    print '\n\n'
  print
  print
  print '''Basis sets implemented in G09.\nSource: http://www.gaussian.com/g_tech/g_ur/m_basis_sets.htm'''
  print
  for key, value in basis.iteritems():
    print '%s : \n' % key
    for basis_set in value:
      print basis_set,
    print '\n\n'
  print 'Diffuse and polarized functions can be added.'
  print