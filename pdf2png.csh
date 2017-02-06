#!/bin/tcsh

if ( $#argv > 0 ) then

  set n = 1
  while ( $n <= $#argv )
    
    inkscape -f $argv[$n] -e ${argv[$n]:r}.png
    @ n ++

  end

else

  set pdfs = `\ls *.pdf`
  
  foreach pdf ( $pdfs )
  
    inkscape -f $pdf -e ${pdf:r}.png
  
  end

endif
