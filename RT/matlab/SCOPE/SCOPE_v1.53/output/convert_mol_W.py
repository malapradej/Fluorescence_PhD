#!/usr/bin/python

'''A script to convert all numbers in a file from umol m-2 s-1 to W m-2.
Arguments: <input_file_name> <header_lines>
<delimiter>
Result: Converts the input file (typically layer_aPAR_Fs_spectral.dat)and 
saves result to output file name with same file name ending in *_wm2.dat. 
2nd argument is the number of header lines, if none use 0. 3rd argument
is the optional delimiters, only use if different from space or tab.  
The conversion is based on the examples in 
http://5e.plantphys.net/article.php?ch=t&id=131
'''

import numpy as np
import sys

ifile = sys.argv[1]
ofile = '%s_wm2.dat' %(ifile.split('.')[0])

if len(sys.argv) == 4:
  delim = sys.argv[3]
  nhead = sys.argv[2]
  arr = np.genfromtxt(ifile, delimiter=delim, skip_header=nhead)

if len(sys.argv) == 3:
  nhead = int(sys.argv[2])
  arr = np.genfromtxt(ifile, skip_header=nhead)

if len(sys.argv) == 2:
  arr = np.genfromtxt(ifile)

start = 400.        # the starting wavelength in nm
nl, nwl = arr.shape
wlarr = np.arange(start, start+nwl)
E = 1.988e-16 / wlarr
E = np.tile(E, (nl, 1))

mol = 10.**-6
avo = 6.02e23

# convert micromoles (umol) per m2 per sec to W per m2
arr = np.multiply(arr*mol*avo, E)

np.savetxt(ofile, arr, delimiter='\t', fmt='%.5f') 
print 'Conversion saved in file:', ofile
