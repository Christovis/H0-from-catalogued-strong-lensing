# Create input files for Gravlens

# Dataset options
los=no_los
nimgs=4
version=b  # version(nimgs=2)=[a,b] ; version(nimgs=2)=[a,b,c]
infile=../limg_catalogs_${los}_${nimgs}_${version}.json
inbase=./Quasars/input/${los}/nimgs_${nimgs} # Location of inputs for gravlens
outbase=./Quasars/results/${los}/nimgs_${nimgs}

# Optimisation parameters
restarts_1=1  # at fixed shear, reoptimize galaxy mass and e/PA
restarts_2=1  # optimize shear along with galaxy mass and e/PA
restarts_3=0  # optimize density slope while keeping rest fixed
restarts_4=3  # optimize everything
restarts_5=1  # optimize H0

# Measurement parameters
pos_error=0.001  #[arcsec]
mu_error=0.05    #[magnitude]
dt_error=0.05    #[hours]

module unload python
module load python/3.6.5

mpirun -np 15 python3 ./quasar_input.py $infile $inbase $outbase $los $restarts_1 $restarts_2 $restarts_3 $restarts_4 $restarts_5 $pos_error $mu_error $dt_error
