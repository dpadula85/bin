#!/bin/tcsh

if ( $#argv > 0 ) then

  set n = 1
  while ( $n <= $#argv )
    
    inkscape -f $argv[$n] -E ${argv[$n]:r}.eps
    @ n ++

  end

else

  set svgs = `\ls *.svg`
  
  foreach svg ( $svgs )
  
    inkscape -f $svg -E ${svg:r}.eps
  
  end

endif
