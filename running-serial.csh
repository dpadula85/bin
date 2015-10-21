#!/bin/tcsh
#
# This script submit eet calculation for gdvh23dev version.
# Place the script where you want and add it to your path.
#
# Execute the script in the folder containing the input files.
# Input files must be organized in the following way:
#
#     1/1.com   \_ site energies
#     2/2.com   /
#     ....
#     V_1.2/V_1.2.com \_ couplings
#     V_1.3/V_1.3.com /
#     ...
#
# In the folder must be present a file called calclist.in containing
# the calculation options. The structure of the file is:
#
# 1st row    >>>   FL_prop FL_coup FL_vacuo FL_chstat
#
#            FL_prop = 0 ... skip site energies calculation
#                    = 1 ... perform site energies calculation
#                                                                         
#            FL_coup = 0 ... skip coupling calculation
#                    = 1 ... perform coupling calculation
#                                                                         
#            FL_vacuo = 0 .. do not perform vacuum calculation
#                     = 1 .. in vacuo calculation wll be performed before of the 
#                            MMPol calculation which used the chk file of vacuo
#                            for initial guess. If FL_vacuo = 1 in the folders 
#                            containing site energies and couplings must be present
#                            the file X_vac.com ( V_X.Y_vac.com)
#                                                                         
#            FL_chstat = 0 . coupling will be computed between all transitions.
#                      = 1 . coupling will be computed only between the most 
#                            intense transitions (file containig the electronic
#                            densities are manipulated by seltran.csh script).
#
#
# next rows  >>>   RES: XXX
#
#                  XXX correspond to the file name for site energy calculation
#                  Note that if you want to exclude some coupling calculations you have
#                  to delete the corresponding folder (if the coupling folder is not 
#                  found the coupling is automatically skipped)
#
# To submit the script on bright cluster type:
#
# qsub -q [QUE] -l nodes=1:ppn=[NPROC]:[RES] -N [JOBNAME] running-serial.csh 
#
#    QUE .....  que name
#    NPROC ...  number of processor that will be used in the local node
#    RES .....  machine type (if que is xxx_old RES should be set = oldres,
#               if que is xxx_new RES should be set = newres)
#    JOBNAME .  name of the job (optional)
#
#  ... enjoy!
#

#
# Folder with the input files:
#
 setenv STDir $PBS_O_WORKDIR
 #setenv Log   $STDir/logbook.$PBS_JOBID
 setenv Log   $STDir/logbook

#echo "EET Calculation"               >> $Log
#echo "Calulations will run on node:" >> $Log
#cat $PBS_NODEFILE                    >> $Log

 if ( ! -e /local/$USER ) mkdir /local/$USER

#
# Set Gaussian environment
#
 setenv vs gdvH23dev
 if ( -e /local/$USER/gauscr.$PBS_JOBID ) rm -r /local/$USER/gauscr.$PBS_JOBID
 mkdir /local/$USER/gauscr.$PBS_JOBID
 setenv gdvroot /cm/shared/apps/gaussian/$vs
 setenv GAUSS_EXEDIR /cm/shared/apps/gaussian/$vs/gdv/
 setenv GAUSS_SCRDIR /local/$USER/gauscr.$PBS_JOBID
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

#echo "Gaussian version: `which gdv`"                                         >> $Log
#echo "prop=$config[1] coup=$config[2] vacuo=$config[3] seltran=$config[4]"   >> $Log
#echo                                                                         >> $Log
#echo "Number of chromophores  : $nres "                                      >> $Log
#echo "List of chromophores    : $reslist "                                   >> $Log
#echo                                                                         >> $Log

#
# Execute Site Energies Calculation on the local machine
#
  if ( $FL_prop == 1 ) then
    foreach i ( $reslist )
      set dgrep1 = "`printf 'PROP %-11s submitted\n' $i`"
      set dgrep2 = "`printf 'PROP %-11s ended successfully\n' $i`"
      set subi1  = `grep "$dgrep1" $Log | wc -l | awk '{print $1}'` 
      set subi2  = `grep "$dgrep2" $Log | wc -l | awk '{print $1}'` 
      if ( $subi1 == 0 && $subi2 == 0 ) then
        set sdate = `date +%d"-"%b"-"%Y"-"%H:%M:%S`
        set sdats = `date +%s`
        printf "PROP %-11s submitted          >       %20s\n" $i $sdate >> $Log
        set dprev = "`printf 'PROP %-11s submitted          >       %20s\n' $i $sdate`"
        cp -r $STDir/$i .
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
        set edate = `date +%d"-"%b"-"%Y"-"%H:%M:%S`
        set edats = `date +%s`
        set eltime = `date -u -d "0 $edats seconds - $sdats seconds" +"%H:%M:%S"`
        set check = `grep "Normal termination" ${i}.log | wc -l | awk '{print $1}'`
        if ( $check == "0" ) then
          set dinfo = "`printf 'PROP %-11s ERROR TERMINATION  > start %20s > end %20s > ((!)) - ((!)) - ((!)) - ((!)) \n' $i $sdate $edate`"
        else
          set dinfo = "`printf 'PROP %-11s ended successfully > start %20s > end %20s > elapsed time %8s  \n' $i $sdate $edate $eltime`"
        endif 
        sed -i "s/$dprev/$dinfo/" $Log
        rm fort.7 CisAttA CisAttB CisDetA CisDetB MMPCMMX MMPMTRX QMTRX QMTRXINF QMTRXT MMPCMMXINF MMPMTRXINF
        cd ..
        cp -r $i $STDir
      endif
    end  
  endif

#
# 4. Execute Couplings Calculations on the local machine
#
  if ( $FL_coup == 1 ) then 
    foreach i ( $reslist )
      set dgrepi = "`printf 'PROP %-11s ended successfully\n' $i`"
      set donei  = `grep "$dgrepi" $Log | wc -l | awk '{print $1}'`
      foreach j ( $reslist )
        set dgrepj = "`printf 'PROP %-11s ended successfully\n' $j`"
        set donej  = `grep "$dgrepj" $Log | wc -l | awk '{print $1}'`
        if ( $j > $i ) then
          if ( $donei == "1" && $donej == "1" ) then
            set dgrepij1 = "`printf 'COUP %-5s %-5s submitted\n' $i $j`"
            set dgrepij2 = "`printf 'COUP %-5s %-5s ended successfully\n' $i $j`"
            set subij1  = `grep "$dgrepij1" $Log | wc -l | awk '{print $1}'`
            set subij2  = `grep "$dgrepij2" $Log | wc -l | awk '{print $1}'`
            if ( -e $STDir/V_$i.$j && $subij1 == 0 && $subij2 == 0 ) then
              set sdate = `date +%d"-"%b"-"%Y"-"%H:%M:%S`
              set sdats = `date +%s`
              printf "COUP %-5s %-5s submitted          >       %20s\n" $i $j $sdate >> $Log
              set dprev = "`printf 'COUP %-5s %-5s submitted          >       %20s\n' $i $j $sdate`"
              cp -r $STDir/V_$i.$j .
              cd V_$i.$j
              if ($FL_chstat == 0 && $FL_vacuo == 0 ) then
                foreach x ( CisPA CisPB CisPATr CisPBTr Dat )
                  cp $STDir/$i/$x ${x:r}
                  cp $STDir/$j/$x ${x:r}2
                end
                gdv V_$i.$j.com

              else if ($FL_chstat == 1 && $FL_vacuo == 0 ) then
                foreach x ( CisPA CisPB CisPATr.sel CisPBTr.sel Dat.sel )
                  cp $STDir/$i/$x ${x:r}
                  cp $STDir/$j/$x ${x:r}2
                end
                gdv V_$i.$j.com

              else if ($FL_chstat == 0 && $FL_vacuo == 1 ) then
                foreach x ( CisPA CisPB CisPATr CisPBTr Dat )
                  cp $STDir/$i/${x}_vac ${x:r}
                  cp $STDir/$j/${x}_vac ${x:r}2
                end
                gdv V_$i.${j}_vac.com
                foreach x ( CisPA CisPB CisPATr CisPBTr Dat )
                  cp $STDir/$i/$x ${x:r}
                  cp $STDir/$j/$x ${x:r}2
                end
                gdv V_$i.$j.com

              else if ($FL_chstat == 1 && $FL_vacuo == 1 ) then
                foreach x ( CisPA CisPB ) 
                  cp $STDir/$i/${x}_vac ${x:r}
                  cp $STDir/$j/${x}_vac ${x:r}2
                end
                foreach x ( CisPBTr CisPATr CisPBTr Dat )
                  cp $STDir/$i/${x}_vac.sel ${x:r}
                  cp $STDir/$j/${x}_vac.sel ${x:r}2
                end 
                gdv V_$i.${j}_vac.com

                foreach x ( CisPA CisPB CisPATr.sel CisPBTr.sel Dat.sel )
                  cp $STDir/$i/$x ${x:r}
                  cp $STDir/$j/$x ${x:r}2
                end
                gdv V_$i.${j}.com
              endif

              set edate = `date +%d"-"%b"-"%Y"-"%H:%M:%S`
              set edats = `date +%s`
              set eltime = `date -u -d "0 $edats seconds - $sdats seconds" +"%H:%M:%S"`
              set check = `grep "Normal termination" V_$i.$j.log | wc -l | awk '{print $1}'`
              if ( $check == "1" ) then
                set dinfo = "`printf 'COUP %-5s %-5s ended successfully > start %20s > end %20s > elapsed time %8s  \n' $i $j $sdate $edate $eltime`"
                sed -i "s/$dprev/$dinfo/" $Log
              endif
              rm Cis* Dat* fort.7 *.chk >& /dev/null
              cd ..
              cp -r V_$i.$j $STDir

            else
              #echo "COUP $i.$j    skipped (FOLDER NOT FOUND) at `date` --- `date +%s`" >> $Log
            endif
          else
            set dinfo = "`printf 'COUP %-5s %-5s skipped            > \n' $i $j `"
          endif
        endif
      end
    end
  endif
  cd .. 
#
# 5. End of script and final cleaning
#
# echo " ................... : DONE      >>> `date +%s`  `date`" >> $Log
# echo >> $Log
  rm -r /local/$USER/gauscr.$PBS_JOBID
  #tcsh cleaning.csh $snap >& /dev/null
#
# 6. Copy back the results
#
#  cd $WDir
#  cp -r * $STDir
#
# END
#
