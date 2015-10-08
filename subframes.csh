#!/bin/tcsh

# Take start and end as parameters
# Modify here their values if preferred
#
# set start = $1
# set end = $2
#
# List of frames
#
set start = 00001 
set end = 00270
set framelist = `seq -w $start $end`

#
# List of residues and couplings in each frame
#
set WDir = `pwd`
set reslist = `awk '{print $2}' $Wdir/reslist.in`
set ResIDi = `awk '{print $1}' $Wdir/couplist.in`
set ResIDj = `awk '{print $2}' $Wdir/couplist.in`

#
# Cycle over frames
#
set init = `expr $start - 1`
foreach frame ( $framelist )

  cd $Wdir/$frame/

  set fin = $frame
  set lim = `expr $frame - $init`

  #
  # Group the frames 4 by 4 and set the correct values of
  # processors, memory and queue for the current frame
  #
  if ( $lim <= 4 ) then
    set nproc = 6
    set mem = 1200
    set q = "short_old"
    set str = $nproc":oldres"
  
  else if ( $lim > 4 && $lim <= 8 ) then
    set nproc = 6
    set mem = 1200
    set q = "long_old"
    set str = $nproc":oldres"
  
  else if ( $lim > 8 && $lim <= 12 ) then
    set nproc = 4
    set mem = 1500
    set q = "short_new"
    set str = $nproc":newres"
  
  else if ( $lim > 12 && $lim < 16 ) then
    set nproc = 4
    set mem = 1500
    set q = "long_new"
    set str = $nproc":newres"

  else if ( $lim == 16 ) then
    set nproc = 4
    set mem = 1500
    set q = "long_new"
    set str = $nproc":newres"
    set init = $frame

  endif
  
  #
  # check from calclist.in which calculations are to be performed
  # and correct residues and/or coupling files according to this
  #
  set config = ( `head -1 calclist.in` )
  set dosite = $config[1]
  set docoup = $config[2]
  
  #
  # nproc and memory correction
  # for each residue in the current frame
  #
  if ( dosite == 1 ) then

    foreach res ( $reslist )

      #
      # Set variables to store line and values of memory and nproc
      #      
      set lineproc = `cat $res/$res.com | grep -i -n nproc | cut -d : -f1`
      set proc = `cat $res/$res.com | grep -i -n nproc | cut -d = -f2`
      set linememfile = `cat $res/$res.com | grep -i -n mem | cut -d : -f1`
      set memfile = `cat $res/$res.com | grep -i -n mem | cut -d = -f2 | sed 's/[A-Za-z]*//g'`

      #
      # Check if the values set in the files correspond to the ones
      # suitable for the queue. If not, change them
      #
      # Check both memory and nproc
      if ( $proc != $nproc  && $memfile != $mem) then
        sed -i "${linemem}s/^%mem=.*/%mem=${mem}MW/" $res/$res.com
        sed -i "${lineproc}s/$proc/$nproc/" $res/$res.com

      # Change processors setting if required
      else if ( $proc != $nproc  && $memfile == $mem) then
        sed -i "${lineproc}s/$proc/$nproc/" $res/$res.com

      # Change memory setting if required
      else if ( $proc == $nproc  && $memfile != $mem) then
        sed -i "${linemem}s/^%mem=.*/%mem=${mem}MW/" $res/$res.com
      endif

    end

  endif

  #
  # nproc and memory correction
  # for each coupling in the current frame
  #
  if ( $docoup == 1 ) then

    set count = 1
    while ( $count <= $#ResIDi )
      
      #
      # Set variables to store line and values of memory and nproc
      #      
      set lineproc = `cat V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com | grep -i -n nproc | cut -d : -f1`
      set proc = `cat V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com | grep -i -n nproc | cut -d = -f2`
      set linememfile = `cat V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com | grep -i -n mem | cut -d : -f1`
      set memfile = `cat V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com | grep -i -n nproc | cut -d = -f2 | sed 's/[A-Za-z]*//g'`

      #
      # Check if the values set in the files correspond to the ones
      # suitable for the queue. If not, change them
      #
      # Check both memory and nproc
      if ( $proc != $nproc  && $memfile != $mem) then
        sed -i "${linemem}s/^%mem=.*/%mem=${mem}MW/" V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com
        sed -i "${lineproc}s/$proc/$nproc/" V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com

      # Change processors setting if required
      else if ( $proc != $nproc  && $memfile == $mem) then
        sed -i "${lineproc}s/$proc/$nproc/" V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com

      # Change memory setting if required
      else if ( $proc == $nproc  && $memfile != $mem) then
        sed -i "${linemem}s/^%mem=.*/%mem=${mem}MW/" V_$ResIDi[$count].$ResIDj[$count]/V_$ResIDi[$count].$ResIDj[$count].com
      endif

      @ count ++
    end

  endif

  #
  # Submit the job with the correct values to the corresponding queue
  #
  qsub -q $q -l nodes=1:ppn=$str -N f$frame ~/bin/running-serial2.csh  

  cd ..

end
