#!/usr/bin/env python
import sys
import numpy as np
#import matplotlib.pyplot as P
import pylab as P
import math
import argparse
from scipy.optimize import leastsq
from scipy.stats import skew,skewtest,mode,shapiro
from scipy.special import erf


# ************************************************************

def basicstat(x):

  xave  = np.average(x)
  sigma = np.std(x)
  nvals = len(x)
  xmin  = min(x) 
  xmax  = max(x) 
  deltax  = xmax - xmin ; spacing = deltax/np.sqrt(nvals)
  nbin = int(deltax/spacing)

  print
  print("N     =  %8d"   % nvals)
  print("Ave   =  %8.4f" % xave)
  print("sigma =  %8.4f" % sigma)
  print("E min =  %8.4f" % xmin)
  print("E max =  %8.4f" % xmax)
  print("N bin =  %3d"   % nbin)
  print

  stat = { 'xave'    : xave ,
           'sigma'   : sigma,
           'xmin'    : xmin ,
           'xmax'    : xmax ,
           'spacing' : spacing,
           'nbin'    : nbin }

  return stat

# ************************************************************

def fitgauskewed(ydata,xdata,average):

  def pdf(x):
    return 1/np.sqrt(2*3.14)*np.exp(-x**2/2)

  def cdf(x):
    return (1 + erf(x/np.sqrt(2))) / 2

  def skewf(p,x):
    e = p[0] 
    w = p[1] 
    a = p[2] 
    t = (x-e) / w
    return 2 / w * pdf(t) * cdf(a*t)

  def errfunc(p,x,y):
    return y - skewf(p,x)

  def maxskew(e,w,a):
    d = a / np.sqrt(1+a**2)
    mean = e + w * d * np.sqrt(2/3.14)
    return mean

  init  = [0.0,1.0,0.0]
  out   = leastsq( errfunc, init, args=(xdata, ydata))
  c = out[0]
  e = c[0] ; w = c[1] ; a = c[2]

  Delta  = max(xdata)-min(xdata)
  XMin   = min(xdata)-(Delta/100)*20
  XMax   = max(xdata)+(Delta/100)*20
  NewX   = np.linspace(XMin,XMax,1000)
  fitted = skewf(c,NewX)

  mean   = maxskew(e,w,a)
  print ("Fitting results:")
  print ("e = %10.4f ; w = %10.4f ; a = %10.4f " % (e,w,a) )
  maximum = max(fitted)
  indmax = np.where(fitted==maximum)
  energymax =  NewX[indmax]
  print ("Maximum of distribution: %8.4f" % energymax)

  return NewX,fitted


def fitgaussian(ydata,xdata,ave,sigma,ngau):


  if   ngau == 1 : 
    fitgau1  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2) 
    errfunc  = lambda p, x, y: (y - fitgau1(p, x))
    init  = [1.0, ave, sigma]

  elif ngau == 2 : 
    fitgau2  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2) + p[3]*np.exp(-0.5*((x-p[4])/p[5])**2)
    errfunc  = lambda p, x, y: (y - fitgau2(p, x))
    init  = [1.0, ave, sigma,  1.0, ave+0.05, 0.02]

  else:
    print("\nNumber of Gaussian out of range\n")
    sys.exit()

  out   = leastsq( errfunc, init, args=(xdata, ydata))

  c = out[0]

  print "Fit Coefficients:"
  print ("Gaussian 1 ... A: %10.4f  mu: %10.4f  sigma: %10.4f " % ( c[0],c[1],abs(c[2])))
  if ngau == 2 :  print ("Gaussian 2.... A: %10.4f  mu: %10.4f  sigma: %10.4f " % ( c[3],c[4],abs(c[5])))
  # wave = (c[1]*c[0]+c[4]*c[3])/(c[0]+c[3])
  # print ("W-Average = %10.7f" % wave)

  Delta = max(xdata)-min(xdata)
  XMin = min(xdata)-(Delta/100)*20
  XMax = max(xdata)+(Delta/100)*20
  NewX  = np.linspace(XMin,XMax,1000)

  if ngau == 1 : fitted = fitgau1(c,NewX)
  if ngau == 2 : fitted = fitgau2(c,NewX)

  return NewX,fitted


# ************************************************************

def hist(x,stat,OPT):

  n, bins, patches = P.hist(x, stat['nbin'],normed=True,stacked=True, histtype='bar')
  P.setp(patches, 'facecolor', 'r', 'alpha', 0.75)

  NData = len(n)


  xint = stat['spacing']
  xnew = np.empty(stat['nbin'])
  for i in range(stat['nbin']):
    xnew[i] = bins[i]+xint/2

  if OPT['fitfun'] == "gau" :
    NewX,fitted = fitgaussian(n,xnew,stat['xave'],stat['sigma'],OPT['ngau'])
  elif OPT['fitfun'] == "skewgau":
    NewX,fitted = fitgauskewed(n,xnew,stat['xave'])
  else:
    print("\nConfused in the choice of fit function\n")
    sys.exit()

#  Dist = np.column_stack((xnew,n))
#  Fit  = np.column_stack((xnew,fitted))
#  np.savetxt('dist.dat',Dist)
#  np.savetxt('fit.dat',Fit)

  P.plot(NewX,fitted, 'k-', linewidth=2)
  P.xlabel('E (eV)', fontsize=14)
  P.ylabel('Counts', fontsize=14)
  P.suptitle('Site energy distribution',fontsize=14)

#  fit1gau  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
#  gau1 = fit1gau(c[0:3],xnew)
#  gau2 = fit1gau(c[3:6],xnew)
#  P.plot(xnew,gau1, 'k-', linewidth=1)
#  P.plot(xnew,gau2, 'k-', linewidth=1)


  # Save data
  if OPT['v'] > 0 :
    cbins = []
    Nbins=len(bins)
    for i in range(Nbins-1):
      cbins.append(((bins[i+1]-bins[i])/2)+bins[i])
    cbins=np.array(cbins)
    Out   = np.empty((NData,2))
    Out[:,0]  = cbins ; Out[:,1]  = n
    OutFile = "hist.%d.dat" % OPT['col']
    np.savetxt(OutFile,Out,fmt="%10.4f")
  
    OutFile = "fit.%d.dat" % OPT['col']
    NFit = len(NewX)
    Out = np.empty((NFit,2))
    Out[:,0]  = NewX ; Out[:,1]  = fitted
    np.savetxt(OutFile,Out,fmt="%10.4f")

    

  if OPT['plot'] : P.show()


# ************************************************************

if __name__ == "__main__":

  # Print the welcome message
  print ("\nMoLECoLab Tools\n > distribution.py  - vs. 1.0\n > Created by Sandro Jurinovich\n > sandro.jurinovich@for.unipi.it\n")

  # Set default options:
  OPT = {}

  # Parse the input string
  parser = argparse.ArgumentParser()
  parser.add_argument("-i",help="Input file")
  parser.add_argument("-v",help="Verbosity",action="count")
  parser.add_argument("-e",help="Exclude the first coulumn",action="store_true",default="False")
  parser.add_argument("-c",help="Number of the coloumn to process",type=int,default="1")
  parser.add_argument("-all",help="Number of the coloumn to process",action="store_true",default="False")
  parser.add_argument("-noplot",help="The program does not show the plot",action="store_false",default="True")
  parser.add_argument("-ngau",help="Number of Gaussian used in the fitting (active when -f = gau)",type=int,default="1")
  parser.add_argument("-f",help="Chose fit function",choices=['gau','skewgau'],default="gau")
  args    = parser.parse_args()
  
  OPT['v']       = args.v
  OPT['plot']    = args.noplot
  OPT['input']   = args.i
  OPT['skipfirst'] = args.e
  OPT['fitfun']  = args.f
  OPT['ngau']    = args.ngau 
  OPT['all']     = args.all
  OPT['col']     = args.c
  Col            = args.c-1

  print("Processing %s" % OPT['input'])
  print("Plotting distribution of data in coloumn: %3d" % (Col+1) )

  # Read the input data 
  DATA = np.loadtxt(OPT['input'])

  if Col < 0 : 
    print ("Error: you cannot select -c < 1")
    sys.exit()

  # Keep the colulumn corresponding to the selected residue
  if ( OPT['skipfirst'] == True ): 
    ChromData = DATA[:,Col+1]
  else:
    ChromData = DATA[:,Col]
 
 
  # ---> TO BE IMPROVED ....
  if OPT['all'] == True:
    OPT['plot'] = False
    for i in range(len(DATA)):
      stat = basicstat(DATA[:,i])
      hist(DATA[:,i],stat,OPT)
  
  else:
    # Performs basic statistical analysis
    stat = basicstat(ChromData)

    # Plot the histrogram
    hist(ChromData,stat,OPT)

  
