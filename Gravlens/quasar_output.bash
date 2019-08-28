echo "Creates .json file of Gravlens results"

# Dataset options
los=$1
nimgs=$2
ext_kappa_file=$3
check=$5
inbase=./Quasars/results/${los}/nimgs_${nimgs}/ # Location of gravlens results

module unload python
module load python/3.6.5

python3 ./quasar_output.py $inbase $los $nimgs $ext_kappa_file $4

if [ "$check" == "yes" ]; then 
    python3 ./quasar_output_analyses.py $inbase $los $nimgs $ext_kappa_file $4
fi 
