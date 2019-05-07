# Create input files for Gravlens

# Dataset options
los=no_los
nimgs=2
infile=../limg_catalogs_${los}_${nimgs}.json
inbase=./Quasars/input/${los}/nimgs_${nimgs} # Location of inputs for gravlens
outbase=./Quasars/results/${los}/nimgs_${nimgs}

# Optimisation options
restars_fitnr_1=1
restars_fitnr_2=1
restars_fitnr_3=3
restars_fitnr_4=1

module unload python
module load python/3.6.5

mpirun -np 15 python3 ./quasar_input.py $infile $inbase $outbase $restars_fitnr_1 $restars_fitnr_2 $restars_fitnr_3 $restars_fitnr_4
