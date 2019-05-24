#!/bin/bash -l
#SBATCH --ntasks=8    # nr. of tasks
#SBATCH -t 0-05:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QU4nLOS 
#SBATCH -o QU4nLOS.out
#SBATCH -e QU4nLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Change to the directory where the job was submitted
module purge
module load gnu_comp/7.3.0 openmpi/3.0.1

# Dataset options
los=no_los
nimgs=4
inbase=./Quasars/input/${los}/nimgs_${nimgs} # Location of inputs for glafic
outbase=./Quasars/results/${los}/nimgs_${nimgs} # Location of results for glafic

rm $outbase/fitH0_*

mpirun -np 8 ./parallel_tasks 0 359 "./glafic ${inbase}/optimize_%d.input"
