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
los=with_los
nimgs=2
version=b  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
infile=../limg_catalogs_${los}_${nimgs}_${version}.json

# Lensing Parameters
dt_sigma=2 # 1-sigma uncertainties in the time-delay measurement (in units of days)
image_amps_sigma=0.3 
flux_ratio_error=0.1
astrometry_sigma=0.004 # 1-sigma astrometric uncertainties of the image positions 
                       # (assuming equal precision)

#module purge
module unload python
module load python/3.6.5

mpirun -np 1 python3 ./quasar_run_092.py $infile $nimgs $los $version $dt_sigma $image_amps_sigma $flux_ratio_error $astrometry_sigma
