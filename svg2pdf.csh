#!/bin/tcsh

if ( $#argv > 0 ) then

  set n = 1
  while ( $n <= $#argv )
    
    inkscape -f $argv[$n] -A ${argv[$n]:r}.pdf
    @ n ++

  end

else

  set svgs = `\ls *.svg`
  
  foreach svg ( $svgs )
  
    inkscape -f $svg -A ${svg:r}.pdf
  
  end

endif
