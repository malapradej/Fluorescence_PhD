#!/bin/bash

# a script that creates the leaf_spectral_sif_#_ratio.dat file in the current 
# directory. The notes on 4/11/14 explaines the idea. We are attempting to
# convert from the leaf_fluorescence_absorption_ratio.dat values to a unitless
# quantity of SIF relative to total irradiance on the leaf. The leaf absorption
# (1 - leaf_spectral_albedo.dat) is multiplied by the SIF at a chosen layer. 
# The only parameter to this script is the layer number between 1 and 60.

filename=leaf_spectral_sif_$1_ratio.dat
paste <(awk 'NR>=241 && NR<=451{print 1-$1}' leaf_spectral_albedo.dat) <(awk -v nr="$1" 'FNR==nr{for(i=1; i<=211; i++) print $i}' layer_fluorescence_absorption_ratio.dat) | awk '{printf "%.7f\n", $1*$2}' > $filename
echo wrote file $filename to current directory.
