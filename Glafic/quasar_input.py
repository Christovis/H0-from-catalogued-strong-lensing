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
        args["pos_error"] = sys.argv[5]
        args["mu_error"] = sys.argv[6]
        args["dt_error"] = sys.argv[7]

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
        system_prior = systems_prior[ii]

        # Source observables file
        data_file_name = args["inbase"] + "/source_obs_" + str(ii) + ".dat"
        text_file = open(data_file_name, "w")
        text_file.write("1 " + str(system["nimgs"]) + " 2.0 0.001 \n")
        for jj in range(system["nimgs"]):
            text_file.write(
                "%s %s %s %s %.2f %.2f %s %s \n"
                % (
                    system["ximg"][jj],
                    system["yimg"][jj],
                    system["mu"][jj],
                    args["pos_error"],
                    abs(system["mu"][jj]) * float(args["mu_error"]),
                    system["delay"][jj],
                    args["dt_error"],
                    str(1),
                )
            )
        text_file.close()

        # Prior file (only for lens possible, not for point-source)
        data_file_name = args["inbase"] + "/prior_" + str(ii) + ".dat"
        text_file = open(data_file_name, "w")
        # velocity disperison [km/sec]
        text_file.write("gauss lens %d %d %.2f %.1f \n" % (1, 1, 50.0, 350.0))
        # ellipticity
        text_file.write("gauss lens %d %d %.2f %.1f \n" % (1, 4, 0.05, 0.5))
        # position-angle
        text_file.write("gauss lens %d %d %.2f %.1f \n" % (1, 5, 5.0, 175.0))
        # core-radius
        text_file.write("range lens %d %d %.2f %.1f \n" % (1, 6, 1.0, 5.0))
        # hubble const.
        text_file.write("range hubble %.2f %.2f \n" % (0.5, 0.9))

        # Optimization file
        template_name = args["templates"] + "/template_optimize.input"
        new_file_name = args["inbase"] + "/optimize_" + str(ii) + ".input"
        copyfile(template_name, new_file_name)
        os.system("sed -i '6s@.*@zl " + str(system_prior["zl"]) + "@' " + new_file_name)
        os.system(
            "sed -i '7s@.*@prefix "
            + args["outbase"]
            + "/fitH0_"
            + str(ii)
            + "@' "
            + new_file_name
        )
        # Define lenses and sources
        os.system(
            "sed -i '25s@.*@  lens sie 100.0 0.0 0.0 1.0 1.0 1.0 0.0"
            + "@' "
            + new_file_name
        )
        os.system("sed -i '26s@.*@  point 2.0 0.0 0.0" + "@' " + new_file_name)
        # Define optimization routine
        os.system("sed -i '31s@.*@  1 1 1 1 1 1 0" + "@' " + new_file_name)
        os.system("sed -i '32s@.*@  0 0 0" + "@' " + new_file_name)
        # Define executione commands
        os.system("sed -i '36s@.*@start_command@' " + new_file_name)
        os.system(
            "sed -i '38s@.*@readobs_point "
            + args["inbase"]
            + "/source_obs_"
            + str(ii)
            + ".dat"
            + "@' "
            + new_file_name
        )
        os.system(
            "sed -i '39s@.*@parprior "
            + args["inbase"]
            + "/prior_"
            + str(ii)
            + ".dat"
            + "@' "
            + new_file_name
        )
        SIEPOI_name = args["outbase"] + "/SIE_POI_" + str(ii)
        os.system("sed -i '41s@.*@optimize point@' " + new_file_name)
        os.system("sed -i '42s@.*@calcein2 2.0 0.0 0.0 1@' " + new_file_name)
        os.system("sed -i '43s@.*@findimg @' " + new_file_name)
        
        # Not Working
        # subprocess.call(["./lensmodel", new_file_name])


if __name__ == "__main__":
    sl_sys_analysis()
