#!/bin/bash
# a script that creates a canopy_albedo.dat file from irradiance_spectra_intercepted.dat and canopy_spectral_outgoing_rad.dat files.
# Only includes values lower than the 953 record due to irradiance values approaching zero at wavelengths greater than this.
# This is due to almost total atmospheric absorption at wavelenghts in the vacinity of 1350nm due to H2O and CO2.

paste <(awk '{print $1}' irradiance_spectra_intercepted.dat) <(awk '{print $1}' canopy_spectral_outgoing_rad.dat)| awk '{if($1==0.) r=0; else r=$2/$1; print r}' | head -953 > canopy_spectral_albedo.dat
echo canopy_spectral_albedo.dat file created in current directory.
