#!/bin/tcsh

if ( $#argv > 0 ) then

  set n = 1
  while ( $n <= $#argv )
    
    inkscape -f $argv[$n] -d 300 -e ${argv[$n]:r}.png
    @ n ++

  end

else

  set svgs = `\ls *.svg`
  
  foreach svg ( $svgs )
  
    inkscape -f $svg -d 300 -e ${svg:r}.png
  
  end

endif
