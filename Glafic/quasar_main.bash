#!/bin/bash -l
#SBATCH --ntasks=20    # nr. of tasks
#SBATCH -t 1-00:00:00  # Runtime in D-HH:MM:SS
#SBATCH -J QU4nLOS
#SBATCH -o ./logs/QU4nLOS.out
#SBATCH -e ./logs/QU4nLOS.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive

# Dataset options
los=no_los
nimgs=4
version=b  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]

# Optimization parameters
priors=1  #0=no 1=yes
opt_1=0
opt_2=0
opt_3=0
opt_4=1
opt_explore=0

# Create .input #################################################################
./quasar_input.bash $los $nimgs $version $priors $opt_1 $opt_2 $opt_3 $opt_4 $opt_explore &&

# Run Gravlens ##################################################################
# nimgs=2 & a -> 999
# nimgs=4 & a -> 359
# nimgs=4 & b -> 471
#./quasar_run_c7.bash $los $nimgs 471 20 &&

# Create .json ##################################################################
./quasar_output.bash $los $nimgs $version $opt_explore

exit
