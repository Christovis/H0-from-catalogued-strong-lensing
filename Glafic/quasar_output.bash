echo "Creates .json file of Glafic results"

# Dataset options
los=$1
nimgs=$2
version=$3  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
inbase=./Quasars/results/${los}/nimgs_${nimgs}/ # Location of gravlens results
opt_explore=$4

module unload python
module load python/3.6.5

python3 ./quasar_output.py $inbase $los $nimgs $version $opt_explore
