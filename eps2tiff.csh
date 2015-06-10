#!/bin/tcsh

foreach x (*.eps)
    set NAME=`echo $x | cut -d '.' -f1`
    convert -density 300 $NAME.eps -resize 1024x1024 -trim +repage $NAME.tiff
end
