#!/bin/tcsh

#
# Default options for CV
# refs and publist are switches that turn on and off
# the presence of references and publication list in the
# cv, and if cv and publications should be on the same file
#
set WDir = "~/Dropbox/CV/latex" 
set refs = 'yes'
set publist = 'yes'

if ( $1 == '' ) then
    set join = 'yes'
else
    set join = $1
endif

if ( $join != 'yes' && $join != 'no' ) then
    echo "The argument should be yes or no."
    exit 0
endif

#
# tex files for cv, and publication list
#
set cv_file = "Padula_cv.tex"
set pub_file = "Padula_pubs.tex"
set pres_file = "Padula_pres.tex"

#
# Screen print of options
#
echo
echo " Options summary"
echo " -----------------------------------"
echo " > References: $refs" 
echo " > Publication list: $publist" 
echo " > Join CV and Publist files: $join"
echo


cd $WDir

#
# References section
# Comment references lines if not already commented
#

# check if lines are commented
set start = `grep -n "References" $cv_file | cut -d : -f1`
set end = `expr $start + 40`
set check = `sed -n "${start}p" $cv_file`

# if yes
if ( $check =~ "%*") then

    # and if references are wanted, delete the comment character
    if ( $refs == 'yes' ) then
        sed -i "${start},${end}s/^%//" $cv_file
    endif

# if not
else

    # and if references are not wanted, comment the lines
    if ( $refs == 'no' ) then
        sed -i "${start},${end}s/^/%/" $cv_file
    endif

endif

#
# Generate cv pdf file
#
pdflatex -interaction=nonstopmode $cv_file > /dev/null

#
# Publications section
# if together with CV, comment intestation line
#
if ( $publist == 'yes' ) then

    if ( $join == 'yes' ) then
        sed -i 's/^\\makecvtitle/%\\makecvtitle/' $pub_file

    else
        sed -i 's/^%\\makecvtitle/\\makecvtitle/' $pub_file

    endif

    pdflatex -interaction=nonstopmode $pub_file > /dev/null
    pdflatex -interaction=nonstopmode $pres_file > /dev/null
    pdfunite ${pub_file:r}.pdf ${pres_file:r}.pdf tmp.pdf

    if ( $join == 'yes' ) then
        pdftk A="${cv_file:r}.pdf" B=tmp.pdf cat A1-1 A2-2 B1-1 B2-2 B3-3 A3-3 output final.pdf
        rm tmp.pdf
        mv final.pdf ${cv_file:r}.pdf
        echo "Your output is in $WDir/${cv_file:r}.pdf"

    else
        mv tmp.pdf Padula_pub.pdf
        echo "Your output is in $WDir/${cv_file:r}.pdf and in $WDir/Padula_pub.pdf"

    endif

endif
