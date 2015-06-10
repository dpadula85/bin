#!/bin/tcsh

foreach x (*.pdf)
    set NAME=`echo $x | cut -d '.' -f1`
    convert -density 300 $NAME.pdf -background white -layers merge $NAME.tiff
end
