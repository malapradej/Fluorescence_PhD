#!/bin/bash
# a script that creates a canopy_spectral_sif_ratio2.dat file from irradiance_spectra_intercepted.dat and canopy_spectral_outgoing_SIF.dat files.
# this script is different from the first in that it normalizes or de-unitizes the sif ratio by integrating the irradiance spectral intercepted over these wavelengths. Thus removes its reliance on the incomming irradiance of individual wavelengths. In essence it normalizes by division by the mean irradiance.
# Only includes 211 records due to Fs starting at 640 to 850 nm. 

mean=$(awk 'NR >= 241 && NR <= 451{tot+=$1}END{print tot/211}' irradiance_spectra_intercepted.dat) 
awk -v mean=$mean '{print $1/mean}' canopy_spectral_outgoing_SIF.dat > canopy_spectral_sif_irr_ratio2.dat
echo canopy_spectral_sif_irr_ratio2.dat file created in current directory.
