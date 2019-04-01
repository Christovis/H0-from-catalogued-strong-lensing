from __future__ import division
import os, sys
from shutil import copyfile
from glob import glob
import numpy as np
import json
import subprocess
from mpi4py import MPI
from mpi_errchk import mpi_errchk

import warnings 
warnings.filterwarnings("ignore", category=RuntimeWarning, append=1)
os.system("taskset -p 0xff %d" % os.getpid())


comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()
mpi_warn_on_fork = 0

@mpi_errchk
def sl_sys_analysis():
    # Get command line arguments
    args = {}
    if comm_rank == 0:
        print(':Registered %d processes' % comm_size)
        args["infile"]  = sys.argv[1]
        args["inbase"]  = sys.argv[2]
        args["outbase"] = sys.argv[3]
        
        # Remove previous input files
        os.system("rm "+args["inbase"]+"/Proc*")
    
    args = comm.bcast(args)

    # Organize devision of strong lensing systems
    with open(args["infile"], 'r') as myfile:
        limg_data = myfile.read()
    systems = json.loads(limg_data)
    sys_nr_per_proc = int(len(systems)/comm_size)
    start_sys = sys_nr_per_proc*comm_rank
    end_sys = sys_nr_per_proc*(comm_rank+1)
    
    if comm_rank == 0:
        print("Each process will have %d systems" % sys_nr_per_proc)
        print("That should take app. %f min." % (sys_nr_per_proc*20))

    for ii in range(len(systems))[start_sys:end_sys]:
        system = systems[ii]
        
        # Write data-file
        data_file_name = args["inbase"] + \
                         "/datafile_"+str(ii)+".dat"
        text_file = open(data_file_name, "w")
        text_file.write("1                # Nr. of lens galaxy \n")
        text_file.write("0.0 0.0 3.0e-03  # position \n")
        text_file.write("0.0 10000.0      # R_eff \n")
        text_file.write("0.0 10000.0      # PA [rad] \n")
        text_file.write("0.0 10000.0      # ellipticity \n")
        text_file.write("1 # nr. of sources \n")
        text_file.write("%s # nr. of source images \n" % system["nimgs"])
        for jj in range(system["nimgs"]):
            text_file.write("%s %s %s %s %s %s %s \n" %
                    (system["ximg"][jj], system["yimg"][jj],
                    system["mags"][jj], str(1e-3), str(1e-3),
                    system["delay"][jj], str(1e-3),
                    )
            )
        text_file.close()
        
        # Write optimization-file
        template_name = args["inbase"]+"/template_optimize.in"
        new_file_name = args["inbase"] + \
                        "/optimize_"+str(ii)+".in"
                        #"/Proc"+str(comm_rank)+"_optimize_"+str(ii)+".in"
        copyfile(template_name, new_file_name)
        os.system(
            "sed -i 's@data.*@data "+data_file_name+"@' "+new_file_name
        )
        fit1_name = args["outbase"]+"/fit1_"+str(ii)
        os.system(
            "sed -i '25s@.*@optimize "+fit1_name+"@' "+new_file_name
        )
        os.system(
            "sed -i '28s@.*@setlens "+fit1_name+".start@' "+new_file_name
        )
        fit2_name = args["outbase"]+"/fit2_"+str(ii)
        os.system(
            "sed -i 's@varyone.*@varyone 1 7 -90 90 37 " + \
            fit2_name+"@' "+new_file_name
        )
        os.system(
            "sed -i '34s@.*@setlens "+fit2_name+".start@' "+new_file_name
        )
        SIEPOI_name = args["outbase"]+"/SIE_POI_"+str(ii)
        os.system(
            "sed -i '37s@.*@optimize "+SIEPOI_name+"@' "+new_file_name
        )
        os.system(
            "sed -i '40s@.*@setlens "+SIEPOI_name+".start@' "+new_file_name
            #"sed -i '40s@.*@ @' "+new_file_name
        )
        fitH0_name = args["outbase"]+"/fitH0_"+str(ii)
        os.system(
            "sed -i '41s@.*@varyh 0.6 0.8 51 "+fitH0_name+"@' "+new_file_name
            #"sed -i '41s@.*@ @' "+new_file_name
        )
        
        #Not Working
        #subprocess.call(["./lensmodel", new_file_name])

if __name__ == "__main__":
    sl_sys_analysis()
