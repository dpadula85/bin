#!/bin/tcsh

if ( $#argv > 0 ) then

  set n = 1
  while ( $n <= $#argv )
    
    inkscape $argv[$n] -d 300 --export-type pdf
    @ n ++

  end

else

  set svgs = `\ls *.svg`
  
  foreach svg ( $svgs )
  
    inkscape $svg -d 300 --export-type pdf
  
  end

endif
