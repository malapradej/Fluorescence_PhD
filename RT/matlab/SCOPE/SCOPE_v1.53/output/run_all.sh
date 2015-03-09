#!/bin/bash

# a script that runs all the other scripts for a certain folder.
# the only argument required is the path to the SCOPE run
# folder under the output/ folder.

if [ "$#" -ne 1 ]; then
    echo "One path name parameter required."
    exit 1
fi

cd $1
../mac2unix_all.sh
../convert_mW_W.sh
../canopy_interceptance.py
../canopy_albedo.sh
../canopy_sif_ratio2.sh
../canopy_sif_model.py
echo executed all scripts for directory:
pwd
cd ..
