#!/bin/bash -l
#SBATCH --ntasks=32    # nr. of tasks
#SBATCH -t 1-00:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QU2nLOS
#SBATCH -o ./logs/QU2nLOS.out
#SBATCH -e ./logs/QU2nLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Dataset options
los=no_los
nimgs=2
version=b  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]

# Optimisation parameters                                                        
opt_explore=0 # random parametere exploration
restarts_1=1  # at fixed shear, reoptimize galaxy mass and e/PA                  
restarts_2=1  # optimize shear along with galaxy mass and e/PA                   
restarts_3=0  # optimize density slope while keeping rest fixed                  
restarts_4=2  # optimize everything                                              
restarts_5=1  # optimize H0

# Create .input #################################################################
./quasar_input.bash $los $nimgs $version $restarts_1 $restarts_2 $restarts_3 $restarts_4 $restarts_5 $opt_explore &&

# Run Gravlens ##################################################################
# nimgs=2 & a,b -> 999
# nimgs=4 & a -> 359
# nimgs=4 & b -> 471
./quasar_run_c7.bash $los $nimgs 1019 32 &&

# Create .json ##################################################################
./quasar_output.bash $los $nimgs $version

exit
