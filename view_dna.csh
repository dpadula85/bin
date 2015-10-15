#!/bin/tcsh

set vmdfile = ~/bin/dna.vmd

set new = $1
set old = `grep -i "mol addfile" $vmdfile | awk '{print $3}'`
set oldline = `grep -in "mol addfile" $vmdfile | cut -d : -f1`

sed -i "${oldline}s/${old}/${new}/" $vmdfile
vmd -e $vmdfile
