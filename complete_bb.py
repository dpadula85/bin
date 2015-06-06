#!/usr/bin/env python

#This Script adds the side chains of aminoacids to the protein analysis obtained with proteincut.py
#for EXAT calculations.
#In a dictionary, one defines the correspondence between the originary residues and their
#side chains.
#The script then looks for all the couplings in which the originary residues are involved
#and creates similar couplings in which the originary residues have been substituted by
#their side chain.
#Finally, couplings between the originary residues and their own side chain are generated.

#import standard Python modules
import os, re, shutil

#import personal modules
import G09_files

#Correspondence between bb index and side chain index
trps = { '0029' : '0130',
         '0064' : '0131',
         '0065' : '0132',
         '0109' : '0133',
         '0112' : '0134'}

#Directories whose name contains "V_" followed by some digits
couplings = sorted(filter(lambda x: re.match('V_\d+', x), os.listdir(os.getcwd())))

#Go through the couplings
for coupling in couplings:
  
  #Couplings are named V_chrom1.chrom2
  chrom1 = re.split('[_.]', coupling)[1]
  chrom2 = re.split('[_.]', coupling)[2]
  
  #if any of the bb residues in the dictionary is involved, check which one and create a new folder
  #substituting the bb residue structure with the corresponding side chain.
  if any([chrom1 in trps.keys(), chrom2 in trps.keys()]):
    
    #Get the options from the coupling file to create a new one with same options
    coup = G09_files.input_file(os.path.join(coupling, '%s.com' % coupling))
    G09_opts = coup.options
    
    if chrom1 in trps.keys() and trps[chrom1] != chrom2:
      
      chrom1 = trps[chrom1]     

      #f1 is associated with chrom2 because chrom1 is now one of the artificially added
      #residues, thus having higher index than chrom2.
      f1 = G09_files.input_file(os.path.join(os.getcwd(), chrom2, '%s.com' % chrom2))         
      f2 = G09_files.input_file(os.path.join(os.getcwd(), chrom1, '%s.com' % chrom1))
      G09_opts['structure'] = f1.structure + f2.structure
      G09_opts['chk'] = 'V_%s.%s' % (f1.name.split('.')[0], f2.name.split('.')[0])
      G09_opts['title'] = 'coupling %s - %s' % (f1.name.split('.')[0], f2.name.split('.')[0])      
      f_coup = G09_files.input_file('V_%s.%s.com' % (f1.name.split('.')[0], f2.name.split('.')[0]), G09_opts)
          
      #Create V_chrom1.chrom2 directory and move the .com file in it
      coup_dir = 'V_%s.%s' % (chrom2, chrom1)
      os.makedirs(coup_dir)
      shutil.move(f_coup.name, os.path.join(coup_dir, f_coup.name))
    
    if chrom2 in trps.keys() and trps[chrom2] != chrom1:
      
      chrom2 = trps[chrom2]      
      f1 = G09_files.input_file(os.path.join(os.getcwd(), chrom1, '%s.com' % chrom1))         
      f2 = G09_files.input_file(os.path.join(os.getcwd(), chrom2, '%s.com' % chrom2))
      G09_opts['structure'] = f1.structure + f2.structure
      G09_opts['chk'] = 'V_%s.%s' % (f1.name.split('.')[0], f2.name.split('.')[0])
      G09_opts['title'] = 'coupling %s - %s' % (f1.name.split('.')[0], f2.name.split('.')[0])            
      f_coup = G09_files.input_file('V_%s.%s.com' % (f1.name.split('.')[0], f2.name.split('.')[0]), G09_opts)
          
      #Create V_chrom1.chrom2 directory and move the .com file in it
      coup_dir = 'V_%s.%s' % (chrom1, chrom2)
      os.makedirs(coup_dir)
      shutil.move(f_coup.name, os.path.join(coup_dir, f_coup.name))
      

#Create the couplings between the originary residue and its side chain.
for bb_res, ind_res in trps.iteritems():
  chrom1, chrom2 = bb_res, ind_res
  f1 = G09_files.input_file(os.path.join(os.getcwd(), chrom1, '%s.com' % chrom1))         
  f2 = G09_files.input_file(os.path.join(os.getcwd(), chrom2, '%s.com' % chrom2))
  G09_opts = f1.options
  G09_opts['funct'] = 'cam-b3lyp'
  G09_opts['job'] = 'td nosymm eet=coup IOp(2/12=3)'
  G09_opts['chk'] = 'V_%s.%s' % (f1.name.split('.')[0], f2.name.split('.')[0])
  G09_opts['title'] = 'coupling %s - %s' % (f1.name.split('.')[0], f2.name.split('.')[0])
  G09_opts['structure'] = f1.structure + f2.structure
  G09_opts['add_opts'] = '\n'      
  f_coup = G09_files.input_file('V_%s.%s.com' % (f1.name.split('.')[0], f2.name.split('.')[0]), G09_opts)
  
  #Create V_chrom1.chrom2 directory and move the .com file in it
  coup_dir = 'V_%s.%s' % (chrom1, chrom2)
  os.makedirs(coup_dir)
  shutil.move(f_coup.name, os.path.join(coup_dir, f_coup.name))
