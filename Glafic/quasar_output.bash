# Creates .json file of Glafic results

# Dataset options
los=with_los
nimgs=4
inbase=./Quasars/results/${los}/nimgs_${nimgs}/ # Location of gravlens results

module unload python
module load python/3.6.5

python3 ./quasar_output.py $inbase $los $nimgs
