#!/usr/bin/python
'''A script that calculates the canopy interceptance according to the 
analytical formula proposed by Stenberg (2007) eq 11 to 13 . The only 
argument to this script is the path to the converted to Wm2 
irradiance_spectra_Wm2.dat file. The irradiance values are corrected for canopy
interceptance and a new corrected file is saved as 
irradiance_spectra_intercepted.dat. Adds the canopy interceptance value (i0) to
the pars_and_input.dat file.
Argument: <path to output>
Creates: <irradiance_spectra_intercepted.dat> <adds canopy_interceptance to
pars_and_input.dat>
'''
import sys
import numpy as np
import os
import p_theory as pt

# get directory path
if len(sys.argv) == 1:
    root = os.getcwd()+'/'
else:
    root = sys.argv[1]
    if root[-1] != '/':
        root = root+'/'

# read LAI value, nr layers and lad parameters
pfilename = 'Parameters/inputdata.txt'
pfile = open(root+pfilename, 'r')
plines = pfile.read()
pfile.close()
i = plines.index('LAI')
LAI = np.float(plines[i+3:i+7])
nr_l = 60
print 'read LAI=%.1f from %s file.' % (LAI, pfilename)

# calculate canopy interceptance (i0)
i0 = pt.interceptance(LAI, 'l', 'pois', (-0.35, -0.15))

# appends i0 to pars_and_input.dat file
pfilename = 'pars_and_input.dat'
pfile = open(root+pfilename, 'r')
plines = pfile.readlines()
pfile.close()
if np.shape(plines)[0] == 1:
    print 'file %s not in unix format.' % pfilename
    raise IOError
plines[0] = plines[0].strip('\n')+'i0\n'
plines[1] = plines[1].strip('\n')+'%12.6f\n' % i0
pfile = open(root+pfilename, 'w')
pfile.writelines(plines)
pfile.close()
print 'wrote i0=%.3f to %s file.' % (i0, pfilename)

# correct the irradiance_spectra_Wm2.dat file and save
pfilename = 'irradiance_spectra_Wm2.dat'
irrad = np.genfromtxt(root+pfilename)
irrad *= i0
pfilename = 'irradiance_spectra_intercepted.dat'
np.savetxt(root+pfilename, irrad, fmt='%.5f')
print 'saved intercepted irrandiance to %s file.' % pfilename