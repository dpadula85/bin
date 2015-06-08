#!/bin/tcsh
#
#
# The FL_ options are read from the calclist.nodeXX file which is generated 
# by distribute.csh script. The values for the four available options have to be written
# in the following order: FL_prop FL_coup FL_vacuo FL_chstat (eg. 1 1 0 1 ).
# Gaussian version should be set in the script (default: gdvH23eet).
#
# OPTIONS:  FL_prop = 0 ... site energies calculation skipped.
#                   = 1 ... site energies calculation will be performed.
#
#           FL_coup = 0 ... coupling calculation skipped. 
#                   = 1 ... coupling calculation will be performed.
#
#           FL_vacuo = 0 .. in vacuo calculation will not be performed.
#                    = 1 .. in vacuo calculation wll be performed before of the 
#                           MMPol calculation which used the chk file of vacuo
#                           for initial guess.
#
#           FL_chstat = 0 . coupling will be computed between all transitions.
#                     = 1 . coupling will be computed only between the most 
#                           intense transitions (file containig the electronic
#                           densities are manipulated by seltran.csh script).


#
# Folder with the input files:
#
 setenv STDir $PBS_O_WORKDIR
 setenv Log   $STDir/logbook.$PBS_JOBID
 echo "EET Calculation"               >  $Log
 echo "Calulations will run on node:" >> $Log
 cat $PBS_NODEFILE                    >> $Log

#
# Set Gaussian environment
#
 setenv vs gdvH23dev
 if ( -e /local/$USER/gauscr ) rm -r /local/$USER/gauscr
 mkdir -p /local/$USER/gauscr
 setenv gdvroot /cm/shared/apps/gaussian/$vs
 setenv GAUSS_EXEDIR /cm/shared/apps/gaussian/$vs/gdv/
 setenv GAUSS_SCRDIR /local/$USER/gauscr
 source $gdvroot/gdv/bsd/gdv.login

#
# Load the input file .com which have to be processed.
#
set list      = "$STDir/calclist.in" 
set config    = (`head -1 $list`)
set FL_prop   = $config[1]
set FL_coup   = $config[2]
set FL_vacuo  = $config[3]
set FL_chstat = $config[4]
set reslist = (`grep "RES:" $list | awk '{print $2}'`)
set nres = ${#reslist}

#
# Create Local Folder and copy files:
#
set local_folder = "eet.$PBS_JOBID"
set WDir = `expr /local/$USER/$local_folder `
if (! -e $WDir) mkdir $WDir
rm -rf $WDir/*
cd $WDir

echo "Gaussian version: `which gdv`"                                         >> $Log
echo "prop=$config[1] coup=$config[2] vacuo=$config[3] seltran=$config[4]"   >> $Log
echo                                                                         >> $Log
echo "Number of chromophores  : $nres "                                      >> $Log
echo "List of chromophores    : $reslist "                                   >> $Log
echo                                                                         >> $Log

#
# Execute Site Energies Calculation on the local machine
#
  if ( $FL_prop == 1 ) then
    foreach i ( $reslist )
      cp -r $STDir/$i .
      echo "PROP $i submitted at `date` --- `date +%s`" >> $Log
      cd $i
      if ( $FL_vacuo == 1 ) then
        gdv ${i}_vac.com
        foreach n (CisPA CisPB CisPATr CisPBTr Dat )
          mv $n ${n}_vac
        end 
        if ($FL_chstat != 0 ) tcsh ../../seltran.csh ${i}_vac 1
        gdv ${i}.com
        if ($FL_chstat != 0 ) tcsh ../../seltran.csh ${i} 0
      else if ($FL_vacuo == 0) then
        gdv ${i}.com
        if ($FL_chstat != 0 ) tcsh ../../seltran.csh ${i} 0
      endif
      set check = `grep "Normal termination" ${i}.log | wc -l | awk '{print $1}'`
      if ( $check == "0" ) then
        echo "PROP $i Error termination"  >> $Log
      else
        echo "PROP $i Done      at `date` --- `date +%s`" >> $Log
      endif 
      foreach rem (fort.7 CisAttA CisAttB CisDetA CisDetB MMPCMMX MMPTRX)
        if (-e $rem ) rm $rem
      end
      cd ..
    end  
  endif
#
# 4. Execute Couplings Calculations on the local machine
#
  if ( $FL_coup == 1 ) then 
    foreach i ( $reslist )
      set donei = `grep "PROP $i Done" $Log | wc -l | awk '{print $1}'  `
      foreach j ( $reslist )
        set donej = `grep "PROP $j Done" $Log | wc -l | awk '{print $1}'  `
        if ( $j > $i && $donei == "1" && $donej == "1" && -d $STDir/V_$i.$j) then
          echo "COUP $i.$j    submitted at `date` --- `date +%s`" >> $Log
          cp -r $STDir/V_$i.$j .
          cd V_$i.$j

	  if ($FL_chstat == 0 && $FL_vacuo == 0 ) then
            foreach x ( CisPA CisPB CisPATr CisPBTr Dat )
              cp ../$i/$x ${x:r}
              cp ../$j/$x ${x:r}2
            end
            gdv V_$i.$j.com

          else if ($FL_chstat == 1 && $FL_vacuo == 0 ) then
            foreach x ( CisPA CisPB CisPATr.sel CisPBTr.sel Dat.sel )
              cp ../$i/$x ${x:r}
              cp ../$j/$x ${x:r}2
            end
            gdv V_$i.$j.com

          else if ($FL_chstat == 0 && $FL_vacuo == 1 ) then
            foreach x ( CisPA CisPB CisPATr CisPBTr Dat )
              cp ../$i/${x}_vac ${x:r}
              cp ../$j/${x}_vac ${x:r}2
            end
            gdv V_$i.${j}_vac.com
            foreach x ( CisPA CisPB CisPATr CisPBTr Dat )
              cp ../$i/$x ${x:r}
              cp ../$j/$x ${x:r}2
            end
            gdv V_$i.$j.com

          else if ($FL_chstat == 1 && $FL_vacuo == 1 ) then
            foreach x ( CisPA CisPB ) 
              cp ../$i/${x}_vac ${x:r}
              cp ../$j/${x}_vac ${x:r}2
            end
            foreach x ( CisPBTr CisPATr CisPBTr Dat )
              cp ../$i/${x}_vac.sel ${x:r}
              cp ../$j/${x}_vac.sel ${x:r}2
            end 
            gdv V_$i.${j}_vac.com

            foreach x ( CisPA CisPB CisPATr.sel CisPBTr.sel Dat.sel )
              cp ../$i/$x ${x:r}
              cp ../$j/$x ${x:r}2
            end
            gdv V_$i.${j}.com
          endif

          set check = `grep "Normal termination" V_$i.$j.log | wc -l | awk '{print $1}'`
          if ( $check == "1" ) echo "COUP $i.$j    Done      at `date` --- `date +%s`" >> $Log
          rm Cis* Dat* fort.7
          cd ..
        else
         # echo "$snap COUP $i.$j    Skipped" >> $logbook
        endif
      end
    end
  endif
  cd .. 
#
# 5. End of script and final cleaning
#
  echo " ................... : DONE      >>> `date +%s`  `date`" >> $Log
  echo >> $Log
#  tcsh cleaning.csh $snap >& /dev/null

  # Copy back
  cd $WDir
  cp -r * $STDir


