#!/bin/bash -l
#SBATCH --ntasks=64    # nr. of tasks
#SBATCH -t 1-00:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QU2bwLOS
#SBATCH -o ./logs/QU2bwLOS.out
#SBATCH -e ./logs/QU2bwLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Dataset options
los=with_los  #[with_los, no_los] include line-of-sight structure? [yes, no]
nimgs=2
version=b  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
infile=limg_catalogs_${los}_${nimgs}_${version}.json
#infile=extKappa_limg_catalogs_${los}_${nimgs}_new_2.json

# Optimisation parameters                                                        
opt_explore=0 # random parametere exploration
restarts_a=1  # at fixed shear, reoptimize galaxy mass and e/PA                  
restarts_b=1  # optimize position angle with shear angle
restarts_c=2  # optimize ellipticity with shear 
restarts_d=0  # analyse degeneracy between ellipticity with shear 
restarts_e=1  # optimize shear along with galaxy mass and e/PA                   
restarts_f=0  # optimize density slope while keeping rest fixed                  
restarts_g=3  # optimize everything                                              
restarts_h=0  # analyse degeneracy between ellipticity with shear 
restarts_i=1  # analyse uncertainty of H0
ext_kappa_file=yes  #[yes, no] use external kappa if known

# Create .input #################################################################
./quasar_input.bash $los $nimgs $version $restarts_a $restarts_b $restarts_c $restarts_d $restarts_e $restarts_f $restarts_g $restarts_h $restarts_i $opt_explore $ext_kappa_file $infile &&

# Run Gravlens ##################################################################
# nimgs=2 & a,b -> 1019
# nimgs=4 & a -> 359
# nimgs=4 & b -> 471
./quasar_run_c7.bash $los $nimgs 1019 64 &&

# Output results in .json #######################################################
check=no  # analyize optimisation resuts
./quasar_output.bash $infile $los $nimgs $ext_kappa_file $check $restarts_a $restarts_b $restarts_c $restarts_d $restarts_e $restarts_f $restarts_g $restarts_h $restarts_i

exit
