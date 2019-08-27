echo "Create input files for Gravlens"
echo "for systems ${1} producing ${2} multi.-images version ${3}"
echo "from the dataset file ${14}"

# Dataset options
los=$1  # with_los
nimgs=$2   # 4
version=$3  # version(nimgs=4)=[a,b] ; version(nimgs=2)=[a,b,c]
infile="../${14}"
if [ "${13}" == "yes" ]; then
    if [ "$los" == "with_los" ]; then
        extkappa=../extKappa.bin
    else
        extkappa=-
    fi
else
    extkappa=-
fi
inbase=./Quasars/input/${los}/nimgs_${nimgs} # Location of inputs for gravlens
outbase=./Quasars/results/${los}/nimgs_${nimgs}

# Optimisation parameters
restarts_1=$4  # at fixed shear, reoptimize galaxy mass and e/PA
restarts_2=$5  # optimize position angle with shear angle
restarts_3=$6  # optimize ellipticity with shear
restarts_4=$7  # optimize shear along with galaxy mass and e/PA
restarts_5=$8  # optimize density slope while keeping rest fixed
restarts_6=$9  # optimize everything
restarts_7=${10}  # analyse degeneracy between ellipticity with shear
restarts_8=${11}  # analyse uncertainty of H0
opt_explore=${12} # random parametere exploration

# Measurement parameters
pos_error=0.001  #[arcsec]
mu_error=0.05    #[magnitude]
dt_error=0.05    #[hours]

module unload python
module load python/3.6.5

mpirun -np 15 python3 ./quasar_input.py $infile $inbase $outbase $los $restarts_1 $restarts_2 $restarts_3 $restarts_4 $restarts_5 $restarts_6 $restarts_7 $restarts_8 $opt_explore $pos_error $mu_error $dt_error $extkappa ${13}
