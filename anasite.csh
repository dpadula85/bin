#!/bin/tcsh

set snaplist = (`seq -w 00001 00275`)

rm *.dat

foreach snap ( $snaplist )
  cd $snap/
  # First Bright State
  grep "nm " $snap.log | awk 'BEGIN { FS = "=" } ; {if ($2 > 0.10) print $0 }' | head -1 | awk '{ print $5 }' >> ../S1.e.dat

  # Second Bright State
  grep "nm " $snap.log | awk 'BEGIN { FS = "=" } ; {if ($2 > 0.10) print $0 }' | head -2 | tail -1 | awk '{ print $5 }' >> ../S2.e.dat
  cd ../
end

seq 1 275 > num

paste num S1.e.dat > boh
mv boh S1.e.dat
paste num S2.e.dat > boh
mv boh S2.e.dat

distribution.py -i S1.e.dat -c 1 -v -e
mv hist.1.dat hist.S1.dat
mv fit.1.dat fit.S1.dat
distribution.py -i S2.e.dat -c 1 -v -e
mv hist.1.dat hist.S2.dat
mv fit.1.dat fit.S2.dat

rm num
