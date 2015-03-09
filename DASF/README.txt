The following is brief description of the commands used in Will's scripts.

PBS commands:
N           - jobname
l walltime  - the maximum time for job
l vmem      - maximum virtual memory for job
j oe        - err and std out to same output file
o /path/file    - the output log file
m n         - no email sent
t 1-10      - 1 to 10 shells run of this script. PBS_ARRAYID env variable contains number of shell

module load /path/file - the module file containing env variables etc. for shell is load

let 'expression=bla' - maths expression takes place

files=(/path/*) - lists files in path and records to files array. To see array use echo ${files[@]}
n_files=${#files[@]} - number of elements in array files.



