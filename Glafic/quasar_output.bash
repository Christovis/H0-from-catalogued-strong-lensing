# Creates .json file of Glafic results

# Dataset options
los=no_los
nimgs=4
version=c  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
inbase=./Quasars/results/${los}/nimgs_${nimgs}/ # Location of gravlens results

module unload python
module load python/3.6.5

python3 ./quasar_output.py $inbase $los $nimgs $version
