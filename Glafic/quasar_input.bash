#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH -t 04:00:00
#SBATCH -J QUinput 
#SBATCH -o QUinput.out
#SBATCH -e QUinput.err
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --exclusive

# Dataset options
los=no_los
nimgs=4
infile=../limg_catalogs_${los}_${nimgs}.json
inbase=./Quasars/input/${los}/nimgs_${nimgs}
outbase=./Quasars/results/${los}/nimgs_${nimgs}

# Optimisation options
restars_fitnr_1=1
restars_fitnr_2=3
restars_fitnr_3=1
restars_fitnr_4=1

#module purge
#module load python/3.6.5

mpirun -np 15 python3 ./quasar_input.py $infile $inbase $outbase $restars_fitnr_1 $restars_fitnr_2 $restars_fitnr_3 $restars_fitnr_4
