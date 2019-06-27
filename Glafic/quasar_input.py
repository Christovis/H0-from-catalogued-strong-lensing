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
        print("Registered %d processes" % comm_size)
        args["infile"] = sys.argv[1]
        args["inbase"] = sys.argv[2]
        args["outbase"] = sys.argv[3]
        args["templates"] = sys.argv[4]
        args["los"] = sys.argv[5]
        args["pos_error"] = sys.argv[6]
        args["mu_error"] = sys.argv[7]
        args["dt_error"] = sys.argv[8]
        args["priors"] = bool(int(sys.argv[9]))
        args["opt_1"] = int(sys.argv[10])
        args["opt_2"] = int(sys.argv[11])
        args["opt_3"] = int(sys.argv[12])
        args["opt_4"] = int(sys.argv[13])
        args["opt_explore"] = int(sys.argv[14])

        # Remove previous input files
        os.system("rm " + args["inbase"] + "/optimize*")
        os.system("rm " + args["inbase"] + "/prior*")
        os.system("rm " + args["inbase"] + "/source_obs_*")

    args = comm.bcast(args)

    # Organize devision of strong lensing systems
    with open(args["infile"], "r") as myfile:
        limg_data = myfile.read()
    systems = json.loads(limg_data)
    sys_nr_per_proc = int(len(systems) / comm_size)
    start_sys = sys_nr_per_proc * comm_rank
    end_sys = sys_nr_per_proc * (comm_rank + 1)
    with open("../lens_catalogs_sie_only.json", "r") as myfile:
        limg_data = myfile.read()
    systems_prior = json.loads(limg_data)

    if comm_rank == 0:
        print("There are in total %d strong-lensing systems" % len(systems))
        print("Each process will have %d systems" % sys_nr_per_proc)
        print("That should take app. %f min." % (sys_nr_per_proc * 20))

    for ii in range(len(systems))[start_sys:end_sys]:

        system = systems[ii]
        #system_prior = systems_prior[ii]

        # Source observables file ###############################################
        source_obs_file = args["inbase"] + "/source_obs_" + str(ii) + ".dat"
        text_file = open(source_obs_file, "w")
        text_file.write("1 " + str(system["nimgs"]) + " 2.0 0.001 \n")
        for jj in range(system["nimgs"]):
            text_file.write(
                "%s %s %s %s %.2f %.2f %s %s \n"
                % (
                    system["ximg"][jj],
                    system["yimg"][jj],
                    abs(system["mu"][jj]),
                    args["pos_error"],
                    abs(system["mu"][jj]) * float(args["mu_error"]),
                    system["delay"][jj],
                    args["dt_error"],
                    str(0),
                )
            )
        text_file.close()

        # Prior file (only for lens possible, not for point-source) #############
        prior_file = args["inbase"] + "/prior_" + str(ii) + ".dat"
        text_file = open(prior_file, "w")
        if args["priors"] > 0:
            # sie:velocity disperison [km/sec]
            text_file.write("range lens %d %d %.2f %.1f \n" % (1, 1, 100.0, 400.0))
            # sie:ellipticity
            text_file.write("range lens %d %d %.3f %.3f \n" % (1, 4, 0.01, 0.95))
            # sie:position-angle
            text_file.write("range lens %d %d %.2f %.1f \n" % (1, 5, -90.0, 90.0))
            # sie:core-radius
            text_file.write("range lens %d %d %.3f %.3f \n" % (1, 6, 0.001, 2.0))
            if args["los"] == "with_los":
                # pert:shear
                text_file.write("range lens %d %d %.3f %.3f \n" % (2, 4, 0.001, 0.9))
                # pert:position-angle
                text_file.write("range lens %d %d %.2f %.1f \n" % (2, 5, -90.0, 90.0))
        if args["los"] == "with_los":
            # pert: match x-pos to lens
            text_file.write("match lens %d %d %d %d %.1f %.1f \n" % (2,2,1,2,1.0,0.0))
            # pert: match y-pos to lens
            text_file.write("match lens %d %d %d %d %.1f %.1f \n" % (2,3,1,3,1.0,0.0))
        # hubble const.
        text_file.write("range hubble %.2f %.2f \n" % (0.2, 1.2))

        # Write optimization file #####################################################
        optimize_fname = args["inbase"] + "/optimize_" + str(ii) + ".in"
        prefix_file = args["outbase"] + "/fitH0_" + str(ii)
        text_file = open(optimize_fname, "w")
        text_file.write("# Primary parameters\n")
        text_file.write("omega 0.3089\n")
        text_file.write("lambda 0.6911\n")
        text_file.write("weos -1.0\n")
        text_file.write("hubble 0.6774\n")
        text_file.write("zl 0.5\n")
        text_file.write("prefix %s\n" % prefix_file)
        text_file.write("xmin     -5.0\n")
        text_file.write("ymin     -5.0\n")
        text_file.write("xmax      5.0\n")
        text_file.write("ymax      5.0\n")
        text_file.write("pix_poi   0.2\n")
        text_file.write("maxlev    6\n")
        text_file.write(" \n")
        text_file.write("# Secondary parameters\n")
        text_file.write("chi2_splane    1\n")
        text_file.write("chi2_checknimg 1\n")
        text_file.write("chi2_usemag    0\n")
        text_file.write("chi2_restart   -1\n")
        #text_file.write("hvary          1\n")
        #text_file.write("ran_seed       -46158\n")
        text_file.write(" \n")
        text_file.write("# Define lenses and sources\n")
        if args["los"] == "with_los":
            text_file.write("startup 2 0 1\n")
        else: 
            text_file.write("startup 1 0 1\n")
        text_file.write("   lens sie 300.0 0.0 0.0 0.5 1.0 0.1 0.0\n")
        if args["los"] == "with_los":
            text_file.write("   lens pert 2.0 0.0 0.0 3.958e-02 5.182e+01 0.0 0.0\n")
        text_file.write("   point 2.0 0.0 0.0\n")
        text_file.write("end_startup\n")
        text_file.write(" \n")
        text_file.write("# Define optimizations\n")
        text_file.write("start_setopt\n")
        text_file.write("   1 1 1 1 1 1 0\n")
        if args["los"] == "with_los":
            text_file.write("   0 0 0 1 1 0 0\n")
        text_file.write("   0 1 1\n")
        text_file.write("end_setopt\n")
        text_file.write(" \n")
        text_file.write("# Start Runs\n")
        text_file.write("start_command\n")
        text_file.write("readobs_point %s\n" % source_obs_file)
        text_file.write("parprior %s\n" % prior_file)
        if args["opt_1"] > 0:
            # 1st Optimization: Mass
            text_file.write("varyone 1 1 100 400 100\n")
        if args["opt_2"] > 0:
            # 2nd Optimization: e/Pa
            text_file.write("varytwo 1 4 0.0 0.95 50 1 5 -90.0 90.0 50\n")
        if (args["opt_3"] > 0) and (system["nimgs"] > 2):
            # 3rd Optimization: shear PA
            text_file.write("varyone 2 5 -90.0 90.0 50\n")
        if args["opt_4"] > 0:
            # 4th Optimization: lens
            text_file.write("optimize\n")
        if args["opt_explore"] > 0:
            text_file.write("opt_explore 1000 1000.0\n")
        text_file.write("varycosmo hubble 0.2 1.2 100\n")
        text_file.write("calcein 2.0\n")
        text_file.write("findimg\n")
        text_file.write("quit\n")


if __name__ == "__main__":
    sl_sys_analysis()
