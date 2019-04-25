#!/bin/bash -l
#SBATCH --ntasks=32    # nr. of tasks
#SBATCH -t 0-05:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QU4wLOS 
#SBATCH -o QU4wLOS.out
#SBATCH -e QU4wLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Change to the directory where the job was submitted
module purge
module load gnu_comp/7.3.0 openmpi/3.0.1

# Dataset options
los=with_los
nimgs=4
inbase=./Quasars/input/${los}/nimgs_${nimgs}

mpirun -np 32 ./parallel_tasks 0 359 "./glafic ${inbase}/optimize_%d.input"
