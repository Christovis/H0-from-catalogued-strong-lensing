#!/bin/bash -l
#SBATCH --ntasks=64    # nr. of tasks
#SBATCH -t 1-00:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QU2bnLOS
#SBATCH -o ./logs/QU2bnLOS.out
#SBATCH -e ./logs/QU2bnLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Dataset options
los=no_los  #[with_los, no_los] include line-of-sight structure? [yes, no]
nimgs=2
version=b  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
infile=limg_catalogs_${los}_${nimgs}_${version}.json
#infile=extKappa_limg_catalogs_${los}_${nimgs}_new_2.json

# Optimisation parameters                                                        
opt_explore=0 # random parametere exploration
restarts_1=1  # at fixed shear, reoptimize galaxy mass and e/PA                  
restarts_2=1  # optimize position angle with shear angle
restarts_3=2  # optimize ellipticity with shear 
restarts_4=1  # optimize shear along with galaxy mass and e/PA                   
restarts_5=0  # optimize density slope while keeping rest fixed                  
restarts_6=3  # optimize everything                                              
restarts_7=1  # analyse degeneracy between ellipticity with shear 
restarts_8=1  # analyse uncertainty of H0
ext_kappa_file=yes  #[yes, no] use external kappa if known

# Create .input #################################################################
#./quasar_input.bash $los $nimgs $version $restarts_1 $restarts_2 $restarts_3 $restarts_4 $restarts_5 $restarts_6 $restarts_7 $restarts_8 $opt_explore $ext_kappa_file $infile &&

# Run Gravlens ##################################################################
# nimgs=2 & a,b -> 1019
# nimgs=4 & a -> 359
# nimgs=4 & b -> 471
./quasar_run_c7.bash $los $nimgs 300 64 &&

# Create .json ##################################################################
#./quasar_output.bash $los $nimgs $ext_kappa_file $infile

exit
