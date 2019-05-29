#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 04:00:00
#SBATCH -J QUinput 
#SBATCH -o QUinput.out
#SBATCH -e QUinput.err
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --exclusive

# Dataset options
los=with_los
nimgs=4
version=c  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
infile=../limg_catalogs_${los}_${nimgs}_${version}.json
inbase=./Quasars/input/${los}/nimgs_${nimgs}
templates=./Quasars/input
outbase=./Quasars/results/${los}/nimgs_${nimgs}

# Measurement parameters
pos_error=0.001  #[arcsec]
mu_error=0.05    #[magnitude]
dt_error=0.05    #[hours]

# Optimization parameters
priors=1  #0=no 1=yes

module unload python
module load python/3.6.5

mpirun -np 1 python3 ./quasar_input.py $infile $inbase $outbase $templates $los $pos_error $mu_error $dt_error $priors
