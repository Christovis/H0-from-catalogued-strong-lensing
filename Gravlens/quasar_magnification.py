# — 1. modelled positions of the source (ys1, ys2)
# — 2. Einstein Radius of each lens (Re)
# — 3. modelled positions of the lensed images (x1_i, x2_i).
# — 4. modelled magnification of the lensed images (mu_i).
# — 5. modelled axis-ration and position angles of the lenses (ql and pa) 

from __future__ import division
import os, sys, glob
from shutil import copyfile
import numpy as np
import pandas as pd
import json

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, append=1)
os.system("taskset -p 0xff %d" % os.getpid())

args = {} 
args["infile"] = sys.argv[1]
args["inbase"] = sys.argv[2]
args["outdirstr"] = sys.argv[3]
args["los"] = sys.argv[4]
args["nimgs"] = sys.argv[5]
outdir = os.fsencode(args["outdirstr"])


os.system("rm "+args["inbase"]+"/magnification_*")

# Load dataset
with open("../lens_catalogs_sie_only.json", 'r') as myfile:
    limg_data = myfile.read()
systems_prior = json.loads(limg_data)

# Collect file names
files_fitH0 = []
file_names = "Quasars/results/%s/nimgs_%s/fitH0_*.best" % \
             (args["los"], args["nimgs"])
for file in glob.glob(file_names):
    files_fitH0.append(file)
# sort by losID
files_fitH0.sort(key=lambda f: int(f.split('_')[-1].split('.')[0]))
ref_losID = np.asarray([ss["losID"] for ss in systems_prior])

# Run through files
counter = 0
for file in files_fitH0:
    losID = file.split('_')[-1].split('.')[0]
    #index = np.where(5 == ref_losID)[0]
    system_prior = systems_prior[int(losID)]
    #print(counter, losID, system_prior["losID"])

    with open(file, "r") as f:
        lines = f.readlines()

        # Run through lines in file
        for index in range(len(lines)):
            
            if 'SOURCE' in lines[index].split():
                source = lines[index+1].split()
                ys1 = float(source[2])
                ys2 = float(source[3])
            
            elif 'h:' in lines[index].split():
                h = float(lines[index].split()[1])

    ## Write lensmodel input-file
    template_name = "/cosma7/data/dp004/dc-beck3/H0_measure/Gravlens/" + \
                    "Quasars/input/template_magnification.in"
    new_file_name = args["inbase"] + \
                    "/magnification_"+str(counter)+".in"
    copyfile(template_name, new_file_name)
    os.system(
            "sed -i '4s@.*@set hval = "+str(h)+"@' "+new_file_name
    )
    os.system(
            "sed -i '5s@.*@set zlens = " + \
            str(system_prior['zl'])+"@' "+new_file_name
    )
    os.system(
            "sed -i '6s@.*@set zsrc = "+str(2.0)+"@' "+new_file_name
    )
    fitH0_name = args["outdirstr"]+"fitH0_"+losID + \
                 ".start"
    os.system(
            "sed -i '9s@.*@setlens "+fitH0_name+"@' "+new_file_name
    )
    mu_name = args["outdirstr"]+"magnification_"+losID+".dat"
    os.system(
            "sed -i '12s@.*@findimg "+str(ys1)+" "+str(ys2)+" "+mu_name+"@' "+new_file_name
    )
    counter += 1


#comm = MPI.COMM_WORLD
#comm_rank = comm.Get_rank()
#comm_size = comm.Get_size()
#mpi_warn_on_fork = 0
#
#@mpi_errchk
#def sl_sys_magnifications():
#    # Get command line arguments
#    args = {} 
#    if comm_rank == 0:
#        print(':Registered %d processes' % comm_size)
#        args["infile"]  = sys.argv[1]
#        args["inbase"]  = sys.argv[2]
#        args["outbase"] = sys.argv[3]
#        args["outdirstr"] = sys.argv[1]
#        args["los"] = sys.argv[2]
#        args["nimgs"] = sys.argv[3]
#        outdir = os.fsencode(args["outdirstr"])
#
#        # Remove previous input files
#        print('args["inbase"]', args["inbase"])
#        os.system("rm "+args["inbase"]+"/magnification*")
#
#    args = comm.bcast(args)

