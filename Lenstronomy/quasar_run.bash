#!/bin/bash -l
#SBATCH --ntasks=32
#SBATCH -t 0-10:00:00
#SBATCH -J LSn4m 
#SBATCH -o LSn4m.out
#SBATCH -e LSn4m.err
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --exclusive

# Dataset options
los=no_los
nimgs=4
infile=../limg_catalogs_${los}_${nimgs}.json

# Optimisation options
restars_fitnr_1=1
restars_fitnr_2=1
restars_fitnr_3=3
restars_fitnr_4=1

#module purge
#module load python/3.6.5

mpirun -np 32 python3 ./quasar_run.py $infile $nimgs $los
