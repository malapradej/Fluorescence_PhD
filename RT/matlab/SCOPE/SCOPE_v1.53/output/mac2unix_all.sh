#!/bin/bash

# a script that converts all *.dat files in the directory from
# which the script was run and all subdirectories from mac to
# unix format.

find . -name '*.dat' |xargs mac2unix
echo done converting all files to unix format.
