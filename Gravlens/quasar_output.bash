# Creates .json file of Glafic results

# Dataset options
los=no_los
nimgs=4
inbase=./Quasars/results/${los}/nimgs_${nimgs}/ # Location of gravlens results

# Optimisation options
restars_fitnr_1=1
restars_fitnr_2=1
restars_fitnr_3=3
restars_fitnr_4=1

module unload python
module load python/3.6.5

python3 ./quasar_input.py $inbase $los $nimgs
