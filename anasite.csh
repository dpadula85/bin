#!/bin/tcsh

set snaplist = (`seq -w 00001 00275`)

rm *.dat

foreach snap ( $snaplist )
  cd $snap/
  # First Bright State
  grep "nm " $snap.log | awk 'BEGIN { FS = "=" } ; {if ($2 > 0.08) print $0 }' | head -1 | awk '{ print $5 }' >> ../S1.e.dat

  # Dipolar Strength First bright state
  # grep -i -A 11 'Ground to excited state transition electric dipole moments' $snap.log | tail -10 | awk '{if ($5 > 0.8) print $5 }' | head -1 >> ../S1.dip.dat

  # Second Bright State
  grep "nm " $snap.log | awk 'BEGIN { FS = "=" } ; {if ($2 > 0.08) print $0 }' | head -2 | tail -1 | awk '{ print $5 }' >> ../S2.e.dat

  # Second Bright State
  grep "nm " $snap.log | awk 'BEGIN { FS = "=" } ; {if ($2 > 0.08) print $0 }' | head -3 | tail -1 | awk '{ print $5 }' >> ../S3.e.dat

  # Dipolar Strength First bright state
  # grep -i -A 11 'Ground to excited state transition electric dipole moments' $snap.log | tail -10 | awk '{if ($5 > 0.8) print $5 }' | head -2 | tail -1 >> ../S2.dip.dat
  cd ../
end

seq 1 275 > num

paste num S1.e.dat > boh
mv boh S1.e.dat
#paste num S1.dip.dat > boh
#mv boh S1.dip.dat
paste num S2.e.dat > boh
mv boh S2.e.dat
paste num S3.e.dat > boh
mv boh S3.e.dat
#paste num S2.dip.dat > boh
#mv boh S2.dip.dat

distribution.py -i S1.e.dat -c 1 -v -e
mv hist.1.dat hist.S1.dat
mv fit.1.dat fit.S1.dat
distribution.py -i S2.e.dat -c 1 -v -e
mv hist.1.dat hist.S2.dat
mv fit.1.dat fit.S2.dat
distribution.py -i S3.e.dat -c 1 -v -e
mv hist.1.dat hist.S3.dat
mv fit.1.dat fit.S3.dat

rm num
