#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH -t 0-06:00:00
#SBATCH -J LSw4m 
#SBATCH -o LSw4m.out
#SBATCH -e LSw4m.err
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --exclusive

# Dataset options
los=no_los
nimgs=4
infile=../limg_catalogs_${los}_${nimgs}.json

# Lensing Parameters
dt_sigma=2 # 1-sigma uncertainties in the time-delay measurement (in units of days)
image_amps_sigma=0.3 
flux_ratio_error=0.1
astrometry_sigma=0.004 # 1-sigma astrometric uncertainties of the image positions 
                       # (assuming equal precision)

#module purge
module unload python
module load python/3.6.5

mpirun -np 1 python3 ./quasar_run.py $infile $nimgs $los $dt_sigma $image_amps_sigma $flux_ratio_error $astrometry_sigma
