#!/bin/tcsh

# Take start and end as parameters
# Modify here their values if preferred
#
# set start = $1
# set end = $2

set WDir = `pwd`
set framelist = `\ls -dltr $WDir/00* | awk -F "/" '{print $NF}' | sort -n` 
set reslist = `awk '{print $2}' $WDir/reslist.in`
set ResIDi = `awk '{print $1}' $WDir/couplist.in`
set check = `expr $#reslist + $#ResIDi`
set Log = $WDir/framelog

if ( -e $Log ) rm $Log

foreach frame ( $framelist )

  if ( -e $frame/logbook ) then
    set control = `grep -c "successfully" $WDir/$frame/logbook`
  else
    set control = '0'
  endif

  if ( $control == $check) then
    echo "Frame $frame -> Done." >> $Log
  else
    set dif = `expr $check - $control`
    echo "Frame $frame -> Incomplete, $dif missing." >> $Log
  endif

end
