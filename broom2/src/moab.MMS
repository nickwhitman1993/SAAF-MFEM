#!/bin/bash
##### These lines are for Moab
#MSUB -l walltime=2:00:00
#MSUB -q pbatch
#MSUB -l procs=80
#MSUB -o MMS.output
#MSUB -l depend=605159

##### These are shell commands
echo -n 'This machine is ';hostname
echo -n 'My jobid is '; echo $SLURM_JOBID
echo 'My path is:' 
echo $PATH
echo ' '
echo '--------------------- Job Output Below ---------------------'
srun -n80 ./MMS -s 8 -e 6 -i 1000 -probEpsilon 3 --no-DSA -d 'sine' -g 3 -r 3 -t 1e-13
