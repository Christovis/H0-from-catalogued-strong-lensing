echo "Run Glafic over ${3} systems"
# Change to the directory where the job was submitted
module purge
module load gnu_comp/7.3.0 openmpi/3.0.1

# Dataset options
los=$1
nimgs=$2
inbase=./Quasars/input/${los}/nimgs_${nimgs} # Location of inputs for gravlens
outbase=./Quasars/results/${los}/nimgs_${nimgs} # Location of results for gravlens

if test -f "$outbase/fitH_0.best"; then
    "removing old results"
    rm $outbase/*
    #rm $outbase/Rein_*
fi

# nimgs=2 & a -> 999 
# nimgs=4 & a -> 359
# nimgs=4 & b -> 471
mpirun -np $4 ./parallel_tasks 0 $3 "./lensmodel ${inbase}/optimize_%d.in" &&
rm $outbase/fit1_* &&
rm $outbase/fit2_* &&
rm $outbase/fit3_*
