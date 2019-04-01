#!/bin/bash -l
# Bash script to execute child subprocesses for each core as defined
# in quasar_run_proc.bash
#
# This file will let a core run over all strong lensing systems as
# they were defined in .Quasars/quasar_input.py

# Properties from parent script
processor=$1
property=$2
previous_results=$3
echo "Property is: ${property} Proc. Nr. is: ${processor}"

# Properties from this script
inbase=./Quasars/input/${property}
outbase=./Quasars/results/${property}
proc_file=${inbase}/Proc${processor}_optimize_*

# Run through files
for filename in $proc_file; do
    # file-nr. corresponds to system-nr.
    filenumber="$(echo $filename | cut -d'_' -f 4 | cut -d'.' -f 1)"

    # If H0 fitted already skip, else do
    if [ -f ${outbase}/Proc${processor}_fitH0_${filenumber}.best ]; then
        continue
    else
        
    fi
done

