#!/bin/bash -l
#SBATCH --ntasks=32    # nr. of tasks
#SBATCH -t 1-00:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QUwithLOS 
#SBATCH -o QUwithLOS.out
#SBATCH -e QUwithLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Change to the directory where the job was submitted
module purge
module load gnu_comp/7.3.0 openmpi/3.0.1

mpirun -np 32 ./parallel_tasks 0 999 "./lensmodel ./Quasars/input/with_los/optimize_%d.in"
