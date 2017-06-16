#!/bin/tcsh

if ( $#argv > 0 ) then

  set n = 1
  while ( $n <= $#argv )
    
    inkscape -f $argv[$n] -d 300 -e ${argv[$n]:r}.png
    @ n ++

  end

else

  set pdfs = `\ls *.pdf`
  
  foreach pdf ( $pdfs )
  
    inkscape -f $pdf -d 300 -e ${pdf:r}.png
  
  end

endif
