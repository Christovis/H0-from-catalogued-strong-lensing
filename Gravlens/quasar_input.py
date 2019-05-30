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
        print(":Registered %d processes" % comm_size)
        args["infile"] = sys.argv[1]
        args["inbase"] = sys.argv[2]
        args["outbase"] = sys.argv[3]
        args["los"] = sys.argv[4]
        args["restart_1"] = sys.argv[5]
        args["restart_2"] = sys.argv[6]
        args["restart_3"] = sys.argv[7]
        args["restart_4"] = sys.argv[8]
        args["restart_5"] = sys.argv[9]
        args["pos_error"] = sys.argv[10]
        args["mu_error"] = sys.argv[11]
        args["dt_error"] = sys.argv[12]

        # Remove previous input files
        print('args["inbase"]', args["inbase"])
        os.system("rm " + args["inbase"] + "/optimize*")
        os.system("rm " + args["inbase"] + "/datafile*")

    args = comm.bcast(args)

    # Organize devision of strong lensing systems
    with open(args["infile"], "r") as myfile:
        limg_data = myfile.read()
    systems = json.loads(limg_data)
    sys_nr_per_proc = int(len(systems) / comm_size)
    start_sys = sys_nr_per_proc * comm_rank
    end_sys = sys_nr_per_proc * (comm_rank + 1)

    # with open("../lens_catalogs_sie_only.json", "r") as myfile:
    #    limg_data = myfile.read()
    # systems_prior = json.loads(limg_data)

    if comm_rank == 0:
        print("Each process will have %d systems" % sys_nr_per_proc)
        print("That should take app. %f min." % (sys_nr_per_proc * 20))

    for ii in range(len(systems))[start_sys:end_sys]:

        system = systems[ii]
        # system_prior = systems_prior[system["losID"]]

        # Write data-file #######################################################
        data_fname = args["inbase"] + "/datafile_" + str(ii) + ".dat"
        text_file = open(data_fname, "w")
        text_file.write("1                # Nr. of lens galaxy \n")
        text_file.write("0.0 0.0 3.0e-03  # position \n")
        text_file.write("0.0 10000.0      # R_eff \n")
        text_file.write("0.0 10000.0      # PA [rad] \n")
        text_file.write("0.0 10000.0      # ellipticity \n")
        text_file.write("1 # nr. of sources \n")
        text_file.write("%s # nr. of source images \n" % system["nimgs"])
        for jj in range(system["nimgs"]):
            text_file.write(
                "%s %s %s %s %s %s %s \n"
                % (
                    system["ximg"][jj],
                    system["yimg"][jj],
                    system["mu"][jj],
                    args["pos_error"],
                    abs(system["mu"][jj]) * float(args["mu_error"]),
                    system["delay"][jj],
                    args["dt_error"],
                )
            )
        text_file.close()

        # Write optimization-file ###############################################
        optimize_fname = args["inbase"] + "/optimize_" + str(ii) + ".in"
        text_file = open(optimize_fname, "w")
        text_file.write("# Cosmological parameters\n")
        text_file.write("set omega = 0.3089\n")
        text_file.write("set lambda = 0.6911\n")
        text_file.write("set hvale = 1.0e2\n")
        text_file.write("set zlens = 0.5\n")
        text_file.write("set zsrc = 2.0\n")
        text_file.write("\n")
        text_file.write("# Lensmodel inputs\n")
        text_file.write("set checkparity = 0\n")
        text_file.write("set omitcore 1.0e-6\n")
        text_file.write("set upenalty 1.0e-4\n")
        text_file.write("set gridflag = 0\n")
        text_file.write("set chimode = 0  # source plane\n")
        text_file.write("data %s\n" % data_fname)
        text_file.write("\n")
        optimized_file_nr = 1
        if int(args["restart_1"]) > 0:
            # First Optimization: at fixed shear, reoptimize galaxy mass and e/PA
            text_file.write("# 1st Opt.: add fixed shear, reoptimize galaxy mass and e/PA\n")
            text_file.write("set restart = %s\n" % args["restart_1"])
            text_file.write("setlens 1 1\n")
            text_file.write("alpha 1.0 0.0 0.0 0.1 10.0 0.0 0.0 0.0 0.0 1.0\n")
            text_file.write("1 0 0 1 1 0 0 0 0 0\n")
            fit_name = args["outbase"] + "/fit%d_%d" % (optimized_file_nr, ii)
            text_file.write("optimize %s\n" % fit_name)
            text_file.write("\n")
            optimized_file_nr += 1
        if int(args["restart_2"]) > 0:
            # Second Optimization: optimize shear along with galaxy mass and e/PA
            text_file.write("# 2nd Opt.: optimize shear along with galaxy mass and e/PA\n")
            text_file.write("set restart = %s\n" % args["restart_2"])
            text_file.write("setlens %s.start\n" % fit_name)
            text_file.write("changevary 1\n")
            text_file.write("1 0 0 1 1 1 1 0 0 0\n")
            fit_name = args["outbase"] + "/fit%d_%d" % (optimized_file_nr, ii)
            text_file.write("varyone 1 7 -90 90 37 %s\n" % fit_name)
            text_file.write("\n")
            optimized_file_nr += 1
        if int(args["restart_3"]) > 0:
            # Third Optimization: optimize density slope
            text_file.write("# 3rd Opt.: optimize density slope\n")
            text_file.write("set restart = %s\n" % args["restart_3"])
            text_file.write("setlens %s.start\n" % fit_name)
            text_file.write("changevary 1\n")
            text_file.write("1 0 0 0 0 0 0 0 0 1\n")
            fit_name = args["outbase"] + "/fit%d_%d" % (optimized_file_nr, ii)
            text_file.write("varyone 1 10 0.5 5 37 %s\n" % fit_name)
            text_file.write("\n")
            optimized_file_nr += 1
        if int(args["restart_4"]) > 0:
            # Fourth Optimization: optimize everything
            text_file.write("# 4th Opt.: optimize everything\n")
            text_file.write("set restart = %s\n" % args["restart_4"])
            text_file.write("setlens %s.start\n" % fit_name)
            text_file.write("changevary 1\n")
            text_file.write("1 1 1 1 1 1 1 0 0 0\n")
            fit_name = args["outbase"] + "/fit%d_%d" % (optimized_file_nr, ii)
            text_file.write("optimize %s\n" % fit_name)
            text_file.write("\n")
            optimized_file_nr += 1
        # Fourth Optimization: H0
        text_file.write("# 5th Opt.: optimize H0\n")
        text_file.write("set restart = %s\n" % args["restart_5"])
        text_file.write("setlens %s.start\n" % fit_name)
        fitH0_name = args["outbase"] + "/fitH0_" + str(ii)  # str(system["losID"])
        text_file.write("varyh 0.5 0.9 101 %s\n" % fitH0_name)
        text_file.write("\n")
        Rein_name = args["outbase"] + "/Rein_" + str(ii)  # str(system["losID"])
        text_file.write("calcRein %d %s\n" % (optimized_file_nr, Rein_name))
        text_file.write("1 1 0.5 2 19\n")
        text_file.write("1 5 -90 90 37\n")
        text_file.write("1 7 -90 90 37\n")
        
        text_file.write("quit\n")


if __name__ == "__main__":
    sl_sys_analysis()
