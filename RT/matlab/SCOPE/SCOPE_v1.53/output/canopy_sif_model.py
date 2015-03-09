#!/usr/bin/python
""" A spectral invariant approach to estimating canopy level SIF ratio as portion
of intercepted radiation. This model requires the p_theory.py module to be in
the pythonpath. It requires as a single argument, a the path to the following 
generated files for each parameter in the model:
canopy_spectral_albedo.dat              Wc
leaf_spectral_albedo.dat                wl
leaf_spectral_fluorescence_mat.dat      fsr
canopy_spectral_sif_irr_ratio.dat       SIFr
wl.dat                                  lam
Parameters/inputdata.txt                LAI, LIDFa, LIDFb
The script calculates the canopy_spectral_albedo_est.dat and 
canopy_spectral_sif_irr_ratio_est.dat which are the spectral invariant estimates
which can be compared with similar estimates from SCOPE and it's RT model approach. 
Arguments: <path_to_files> <#: see above>
Returns: <canopy_spectral_albedo_est.dat> <canopy_spectral_sif_irr_ratio_est.dat>
    also plots the graphs of the comparison with SCOPE estimates.
"""

import sys
import os
import p_theory as pt
import numpy as np
import matplotlib.pylab as plt
plt.ion()

fname_alb = 'canopy_spectral_albedo_est.dat'
fname_sif = 'canopy_spectral_sif_irr_ratio_est.dat'

# reads in path
if len(sys.argv)==2:
    root = sys.argv[1]
else:
    root = os.curdir

# read parameters from local files
Wc = np.genfromtxt(os.path.join(root, 'canopy_spectral_albedo.dat'))
wl = np.genfromtxt(os.path.join(root, 'leaf_spectral_albedo.dat'))
fsr = np.genfromtxt(os.path.join(root, 'leaf_spectral_fluorescence_mat.dat'))
SIFr = np.genfromtxt(os.path.join(root, 'canopy_spectral_sif_irr_ratio2.dat'))
lam = np.genfromtxt(os.path.join(root, 'wl.dat'),skip_header=2)

# reads in parameters from the inputdata.txt file
filename = os.path.join(root,'Parameters/inputdata.txt')
f = open(filename, 'r')
lines = f.read()
f.close()
i = lines.index('LAI')
LAI = np.float(lines[i+3:i+7])
i = lines.index('LIDFa')
LIDFa = np.float(lines[i+5:i+11])
i = lines.index('LIDFb')
LIDFb = np.float(lines[i+5:i+11])

# calculate p recollision probability
p = pt.p(LAI, arch='l', par=(LIDFa, LIDFb))

# appends p to pars_and_input.dat file
pfilename = 'pars_and_input.dat'
pfile = open(os.path.join(root,pfilename), 'r')
plines = pfile.readlines()
pfile.close()
if np.shape(plines)[0] == 1:
    print 'file %s not in unix format.' % pfilename
    raise IOError
plines[0] = plines[0].strip('\n')+'\tp\n'
plines[1] = plines[1].strip('\n')+'%12.6f\n' % p
pfile = open(os.path.join(root,pfilename), 'w')
pfile.writelines(plines)
pfile.close()
print 'wrote p=%.6f to %s file.' % (p, pfilename)

# calculate the spectral invariant estimated canopy albedo
We = pt.Alb(wl, p)
np.savetxt(os.path.join(root,fname_alb), We, fmt='%.7f')

# plot comparison between estimated canopy and SCOPE canopy albedo
beg = 0
end = 451
plt.title(r'Canopy and Leaf Spectral Albedo', fontsize=16)
plt.plot(lam[beg:end], wl[beg:end], label='PROSPECT Leaf Albedo')
plt.plot(lam[beg:end], We[beg:end], label='p-theory Canopy Albedo')
plt.plot(lam[beg:end], Wc[beg:end], label='SCOPE Canopy Albedo')
plt.xlabel(r'Wavelength in nm ($\lambda$)', fontsize=14)
plt.ylabel('Albedo', fontsize=14)
plt.legend(loc=2)
plt.savefig(os.path.join(root,'figures','Albedo_leaf_canopy.png'))
plt.clf()

# calculate the spectral invariant estimate of canopy SIF ratio
beg = 240
end = 451
SIFre = pt.SIFr_approx(wl[beg:end], p, fsr)
np.savetxt(os.path.join(root,fname_sif), SIFre, fmt='%.7f')

# plot comparison between estimated canopy and SCOPE canopy SIF ratio
beg = 240
end = 451
plt.title('Canopy and Leaf Spectral SIF ratio', fontsize=16)
plt.plot(lam[beg:end], fsr, label='SCOPE Leaf SIF ratio')
plt.plot(lam[beg:end], SIFre, label='p-theory Canopy SIF ratio')
plt.plot(lam[beg:end], SIFr, label='SCOPE Canopy SIF ratio')
plt.xlabel(r'Wavelength in nm ($\lambda$)', fontsize=14)
plt.ylabel('SIF ratio', fontsize=14)
plt.legend(loc=2)
plt.savefig(os.path.join(root,'figures','SIF_ratio_leaf_canopy.png'))
plt.clf()
plt.close()

print '''created plots in figures/ folder and saved canopy_spectral_albedo_est.dat
and canopy_spectral_sif_irr_ratio_est.dat in current folder.'''