echo "Create input files for Gravlens"
echo "for systems ${1} producing ${2} multi.-images version ${3}"

# Dataset options
los=$1
nimgs=$2
version=$3  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
infile=../limg_catalogs_${los}_${nimgs}_${version}.json
inbase=./Quasars/input/${los}/nimgs_${nimgs}
templates=./Quasars/input
outbase=./Quasars/results/${los}/nimgs_${nimgs}

# Measurement parameters
pos_error=0.001  #[arcsec]
mu_error=0.05    #[magnitude]
dt_error=0.05    #[hours]

# Optimization parameters
priors=$4  #0=no 1=yes
opt_1=$5
opt_2=$6
opt_3=$7
opt_4=$8
opt_explore=$9

module unload python
module load python/3.6.5

mpirun -np 1 python3 ./quasar_input.py $infile $inbase $outbase $templates $los $pos_error $mu_error $dt_error $priors $opt_1 $opt_2 $opt_3 $opt_4 $opt_explore
