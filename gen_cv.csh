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
set join = 'no'

#
# tex files for cv, and publication list
#
set cv_file = "$WDir/Padula_cv.tex"
set pub_file = "$WDir/Padula_pubs.tex"
set pres_file = "$WDir/Padula_pres.tex"

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
set check = `sed -n '241p' $cv_file`

# if yes
if ( $check =~ "%*") then

    # and if references are wanted, delete the comment character
    if ( $refs == 'yes' ) then
        sed -i '241,281s/^%//' $cv_file
    endif

# if not
else

    # and if references are not wanted, comment the lines
    if ( $refs == 'no' ) then
        sed -i '241,281s/^/%/' $cv_file
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
        pdfunite ${cv_file:r}.pdf tmp.pdf final.pdf
        rm tmp.pdf
        mv final.pdf ${cv_file:r}.pdf
        echo "Your output is in ${cv_file:r}.pdf"

    else
        mv tmp.pdf Padula_pub.pdf
        echo "Your output is in ${cv_file:r}.pdf and in Padula_pub.pdf"

    endif

endif
