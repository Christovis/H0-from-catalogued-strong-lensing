#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 04:00:00
#SBATCH -J QUinput 
#SBATCH -o QUinput.out
#SBATCH -e QUinput.err
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --exclusive

property=no_los
inbase=./Quasars/input/${property}
infile=${inbase}/limg_catalogs_${property}_2.json
outbase=./Quasars/results/${property}

#module purge
#module load python/3.6.5

mpirun -np 16 python ./Quasars/quasar_input.py $infile $inbase $outbase
