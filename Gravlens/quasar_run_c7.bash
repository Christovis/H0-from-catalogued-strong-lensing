#!/bin/bash -l
#SBATCH --ntasks=32    # nr. of tasks
#SBATCH -t 1-00:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QU4wLOS 
#SBATCH -o ./logs/QU4wLOS.out
#SBATCH -e ./logs/QU4wLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Change to the directory where the job was submitted
module purge
module load gnu_comp/7.3.0 openmpi/3.0.1

# Dataset options
los=with_los
nimgs=4
inbase=./Quasars/input/${los}/nimgs_${nimgs} # Location of inputs for gravlens
outbase=./Quasars/results/${los}/nimgs_${nimgs} # Location of results for gravlens

rm $outbase/fitH0_*
rm $outbase/Rein_*

# nimgs=2 & a -> 999 
# nimgs=4 & a -> 359
# nimgs=4 & b -> 471
mpirun -np 32 ./parallel_tasks 0 471 "./lensmodel ${inbase}/optimize_%d.in"
rm $outbase/fit1_*
rm $outbase/fit2_*
