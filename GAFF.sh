#!/bin/bash -l

# Instructions to generate parameters force field parameters for 
# a new molecule with "antechamber", a simple routine belonging to 
# the AMBER packages. The automatic parametrisation is based on the 
# Generalized AMBER Force Field technology (GAFF). 

# ****************** #
# ENVIRONMENT SETUP  #
# ****************** #
# All the programs needed are basically AMBER plugins.
# They can be installed by conda, but are already available on my paths.
# Source them using 

source /home/dpadula/Data/software/anaconda3/etc/profile.d/conda.sh
conda activate /home/dpadula/Data/software/anaconda3/envs/AmberTools

# ****************** #
# INPUT              #
# ****************** #
# The input for everything that is coming below, is a SINGLE file. 
# Specifically a Gaussian LOG file where the electrostatic potential (ESP) has been fitted on a 
# grid. In the example below, my input file will be called ---> ${root}.log 
# Please, rename accordingly.
# p.s. ---> ESP and [R]ESP are completely different things. 


# ****************** #
# PARAMETRISATION    #
# ****************** #

# All the commands below are more or less automatic.
root=MAMHPA
antechamber -i ${root}.log -fi gout -gv 0 -o ${root}.mol2 -fo mol2 -c resp -s 2 -at gaff2 -pl 15
parmchk2 -i ${root}.mol2 -f mol2 -o ${root}.frcmod

# -------------------------------------------------------------------------------------- #
#  TLEAP                                                                                 #
#  The part below will be executed within an interactive program (tleap).                #
#                                                                                        #
# The tleap details can be found at: http://ambermd.org/tutorials/basic/tutorial4b/      #
# -------------------------------------------------------------------------------------  #

cat << EOF > ${root}.tleap
source oldff/leaprc.ff94
source leaprc.gaff
MyMOL = loadmol2 ${root}.mol2
loadamberparams ${root}.frcmod
saveamberparm MyMOL ${root}.prmtop ${root}.inpcrd
quit
EOF

tleap -f ${root}.tleap

# The following will generate a folder called "MOL.amb2gmx"

acpype -p ${root}.prmtop -x ${root}.inpcrd

conda activate /home/dpadula/Data/software/anaconda3/envs/working
