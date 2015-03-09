#!/usr/bin/python
"""
A script to create a Fs ratio file which requires 2 matrix files as inputs.
The single argument is the path to the absorption file (commonly 
layer_aPAR_Fs_spectral_wm2.dat), and the file of spectral fluorescence 
(commonly layer_fluorescenceEm_spectral.dat). Both are in units of W m-2 nm-1.
The script will also plot the ratio per layer, per wavelength and save it in the
plots subfolder.
Arguments: <path to files>
Returns: Creates the plots in plots subfolder, and creates a matrix file of Fs 
ratios (commonly layer_Fs_ratios.txt) in same path as the passed files.
"""

import sys
import os
import numpy as np
import matplotlib.pylab as plt

start_wla = 400
end_wla = 850
start_wlf = 640
end_wlf = 850

wla = np.arange(start_wla, end_wla + 1)
wlf = np.arange(start_wlf, end_wlf + 1)
len_a = len(wla)
len_f = len(wlf)

path = sys.argv[1]
if not(os.path.exists(path)):
    print 'path: ', path, 'does not exit.'
    exit()

if path[-1] != '/':
    path = '%s/' %path

abs_file = '%slayer_aPAR_Fs_spectral_wm2.dat' %path
fs_file = '%slayer_fluorescenceEm_spectral.dat' %path
ratio_file = '%slayer_fluorescence_absorption_ratio.dat' %path
plt_file1 = '%sfigures/fluorescence_ratios_plot.png' %path
plt_file2 = '%sfigures/fluorescence_ratios_image.png' %path

abs_arr = np.genfromtxt(abs_file)
fs_arr = np.genfromtxt(fs_file)

# remove soil layer from fluorescence file
if abs_arr.shape[0] < fs_arr.shape[0]:
    fs_arr = fs_arr[:abs_arr.shape[0], :]

# include only common wavelengths in both arrays
skip_wla = start_wlf - start_wla
abs_arr = abs_arr[:, skip_wla:]

# divide arrays to get ratio
ratio_arr = fs_arr / abs_arr
lyr = np.arange(1, ratio_arr.shape[0]+1)

# save to file
np.savetxt(ratio_file, ratio_arr, fmt='%.6f', delimiter='\t')

# plot the ratios to image file
# plot top and bottom canopy ratios per wavelength
fig, ax1 = plt.subplots()
fig.suptitle('Ratio of Fs/absorption, and absorption', fontsize=18)
ax1.plot(wlf, ratio_arr[0], 'b-', label='ratio TOC')
ax1.plot(wlf, ratio_arr[-1], 'b--', label='ratio BOC')
ax1.set_xlabel(r'wavelength in nm $(\lambda)$')
ax1.set_ylabel('Ratio Fs / Absorption', color='b')
ax1.legend(loc=6)

ax2 = ax1.twinx()
ax2.plot(wlf, abs_arr[0], 'g-', label='absorp. TOC')
ax2.plot(wlf, abs_arr[-1], 'g--', label='absorp. BOC')
ax2.set_ylabel(r'$W m^{-2} $', color='g', size='large')
ax2.legend(loc=7)

plt.savefig(plt_file1, format='png')
plt.clf()

# plot the image of the ratios per canopy layer and per wavelength
xx, yy = np.meshgrid(wlf, lyr)
plt.pcolormesh(xx, yy, ratio_arr)
plt.title('Ratio of Fs/absorption', fontsize=14)
plt.axis('image')
plt.gca().invert_yaxis()
plt.ylabel('Layer level')
plt.xlabel(r'wavelength in nm $(\lambda)$')
plt.colorbar(shrink=0.3)
plt.savefig(plt_file2, format='png')
plt.close()

print 'plots saved in figures/ directory, and ratio array saved.'