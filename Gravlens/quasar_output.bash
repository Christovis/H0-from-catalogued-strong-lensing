echo "Creates .json file of Gravlens results"

# Dataset options
infile=$1
los=$2
nimgs=$3
ext_kappa_file=$4
check=$5
inbase=./Quasars/results/${los}/nimgs_${nimgs}/ # Location of gravlens results
restarts_a=$6
restarts_b=$7
restarts_c=$8
restarts_d=$9
restarts_e=${10}
restarts_f=${11}
restarts_g=${12}
restarts_h=${13}
restarts_i=${14}

module unload python
module load python/3.6.5

python3 ./quasar_output.py $inbase $infile $los $nimgs $ext_kappa_file $restarts_a $restarts_b $restarts_c $restarts_d $restarts_e $restarts_f $restarts_g $restarts_h $restarts_i 

if [ "$check" == "yes" ]; then 
    python3 ./quasar_output_analyses.py $inbase $los $nimgs $ext_kappa_file
fi 




