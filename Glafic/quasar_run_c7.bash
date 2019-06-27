echo "Run Glafic over ${3} systems"

# Change to the directory where the job was submitted
module purge
module load gnu_comp/7.3.0 openmpi/3.0.1

# Dataset options
los=$1
nimgs=$2
inbase=./Quasars/input/${los}/nimgs_${nimgs} # Location of inputs for glafic
outbase=./Quasars/results/${los}/nimgs_${nimgs} # Location of results for glafic

rm $outbase/fitH0_*

# nimgs=2 & a -> 999
# nimgs=4 & a -> 359
# nimgs=4 & b -> 471
mpirun -np $4 ./parallel_tasks 0 $3 "./glafic ${inbase}/optimize_%d.in"
