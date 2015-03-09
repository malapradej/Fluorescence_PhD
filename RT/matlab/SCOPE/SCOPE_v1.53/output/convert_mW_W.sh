#!/bin/bash
# a script to convert the irradiance_spectra.dat file from mWm-2 to Wm-2.
# run the script from inside the directory where irradiance_spectra.dat file
# is found.

awk 'FNR==3{for(i=1; i<=NF; i++) print $i/1000}' irradiance_spectra.dat > irradiance_spectra_Wm2.dat
echo saved the irradiance_spectra_Wm2.dat file in current directory.
