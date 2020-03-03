#!/bin/tcsh

if ( $#argv > 0 ) then

  set n = 1
  while ( $n <= $#argv )
    
    inkscape --pdf-poppler $argv[$n] -d 300 --export-type png

    @ n ++

  end

else

  set pdfs = `\ls *.pdf`
  
  foreach pdf ( $pdfs )
  
    inkscape --pdf-poppler $pdf -d 300 --export-type png
  
  end

endif
