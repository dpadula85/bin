#!/usr/bin/env python

#Dictionary for energy conversion
#The key should be of the form u_a2u_b, where u_a and u_b are two different
#units


energy_conversion = {'au'       : 1,
                     'eV'       : 27.21138505,
                     'wn'       : 219474.63,
                     'kj/mol'   : 2625.5,
                     'kcal/mol' : 627.503}
                     
#au to eV = 27.21128505 / 1
#eV to au = 1 / 27.21128505
#u1 to u2 = dict[u2]/dict[u1]

                     
#Dictionaries comparison:
#Options are stored in dictionaries. dictA is the default options dictionary
#and dictB is the one passed to the function.
#We want to compare dictA and dictB. We want to add to dictB all the missing keys
#in dictA with their value.

def dict_compare(dictA, dictB):

  for k in dictA.viewkeys() - dictB.viewkeys():
    dictB[k]= dictA[k]
        
  return dictB                   
