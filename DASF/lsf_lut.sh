#!/bin/bash

# A script that will read in the input lut and loop over every entry and 
# run a bsub command line operation to send jobs to lotus.
# run using: ./lsf_lut.sh <lut_file.dat>
# with <arg>...
WRK_DIR=/work/scratch/malaprade/DASF/Data
COD_DIR=/home/users/malapradej/DASF/Code

DEFAULT_JOB_CWD =$WRK_DIR/LUT

while IFS= read -r line
do
    [[ $line = \#* ]] && continue # skips all # lines
    jbnm=${line// /_}
    echo $fname
    
    bsub -q lotus -n 3 -We 1 -eo $WRK_DIR/err/$jbnm.err -J $jbnm ./Emulator_LUT.py $line
done <$1
