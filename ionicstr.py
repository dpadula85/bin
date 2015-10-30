#!/usr/bin/env python

import argparse as arg

NA = 6.022E23

if __name__ == '__main__':

  parser = arg.ArgumentParser(description="Calculate number of ions to get the ionic strength")
  parser.add_argument('--ionicstrenght','-i',help='Ionic strenght',required=True,type=float)
  parser.add_argument('--boxdimension','-b',help='Box dimension ( l x l x l (Ang))',required=True,nargs=3,type=float)
  args = parser.parse_args()

  IonStr  = args.ionicstrenght
  BoxDim  = args.boxdimension

  BoxVol  = BoxDim[0]*BoxDim[1]*BoxDim[2]

  print " > Ionic strenght (mol/L) : %10.4f" % IonStr
  print " > Box Dimension (Ang)    : %10.4f %10.4f %10.4f" % tuple(BoxDim)
  print " > Box Voulume (Ang)^3    : %10.4f"  % BoxVol

  NMolxL = IonStr*NA     # Particles per liter
  BoxVol = BoxVol*1E-27  # Box Volume in dm^3
  NMolxBox = NMolxL*BoxVol

  print " > Numer of particles     : %d"  % round(NMolxBox)


