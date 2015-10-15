#!/bin/tcsh

set vmdfile = ~/bin/dna.vmd
set prmtop = $1
set netcdf = $2

if ( $netcdf == "" ) then
  echo "Two arguments are required: a prmtop file and a netcdf file"
  exit 0
else
  set oldprmtop = `grep -i "mol new" $vmdfile | awk '{print $3}'`
  set oldprmtopline = `grep -in "mol new" $vmdfile | cut -d : -f1`
  set oldnetcdf = `grep -i "mol addfile" $vmdfile | awk '{print $3}'`
  set oldnetcdfline = `grep -in "mol addfile" $vmdfile | cut -d : -f1`
  
  sed -i -e "${oldprmtopline}s/${oldprmtop}/${prmtop}/ ; ${oldnetcdfline}s/${oldnetcdf}/${netcdf}/" $vmdfile
  vmd -e $vmdfile
endif
