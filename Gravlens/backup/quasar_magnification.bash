# Creates .json file of Glafic results

# Dataset options
los=with_los
nimgs=4
infile=../limg_catalogs_${los}_${nimgs}.json
inbase=./Quasars/input/${los}/nimgs_${nimgs} 
outbase=./Quasars/results/${los}/nimgs_${nimgs}/ # Location of gravlens results

module unload python
module load python/3.6.5

python3 ./quasar_magnification.py $infile $inbase $outbase $los $nimgs
