#!/bin/bash
# a script that creates a canopy_spectral_sif_ratio.dat file from irradiance_spectra_intercepted.dat and canopy_spectral_outgoing_SIF.dat files.
# Only includes 211 records due to Fs starting at 640 to 850 nm. 

paste <(awk 'NR >= 241 && NR <= 451' irradiance_spectra_intercepted.dat) <(awk '{print $1}' canopy_spectral_outgoing_SIF.dat)| awk '{r=$2/$1; print r}' > canopy_spectral_sif_irr_ratio.dat
echo canopy_spectral_sif_irr_ratio.dat file created in current directory.
