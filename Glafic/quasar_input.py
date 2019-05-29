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
        args["priors"] = bool(sys.argv[8])

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
                    abs(system["mu"][jj]),
                    args["pos_error"],
                    abs(system["mu"][jj]) * float(args["mu_error"]),
                    system["delay"][jj],
                    args["dt_error"],
                    str(0),
                )
            )
        text_file.close()

        # Prior file (only for lens possible, not for point-source)
        data_file_name = args["inbase"] + "/prior_" + str(ii) + ".dat"
        text_file = open(data_file_name, "w")
        # sie:velocity disperison [km/sec]
        text_file.write("gauss lens %d %d %.2f %.1f \n" % (1, 1, 300.0, 150.0))
        if args["priors"] is True:
            # sie:ellipticity
            text_file.write("gauss lens %d %d %.2f %.1f \n" % (1, 4, 0.5, 0.4))
            # sie:position-angle
            text_file.write("gauss lens %d %d %.2f %.1f \n" % (1, 5, 175.0, 170.0))
            # sie:core-radius
            text_file.write("range lens %d %d %.2f %.1f \n" % (1, 6, 0.0, 1.0))
        # pert: match x-pos to lens
        text_file.write("match lens %d %d %d %d %.1f %.1f \n" % (2,2,1,2,1.0,0.0))
        # pert: match y-pos to lens
        text_file.write("match lens %d %d %d %d %.1f %.1f \n" % (2,3,1,3,1.0,0.0))
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
        os.system("sed -i '24s@.*@startup 2 0 1" + "@' " + new_file_name)
        os.system(
            "sed -i '25s@.*@  lens sie 300.0 0.0 0.0 0.5 1.0 0.1 0.0"
            + "@' "
            + new_file_name
        )
        os.system(
            "sed -i '26s@.*@  lens pert 2.0 0.0 0.0 3.958251e-02 5.182124e+01 0.0 0.0"
            + "@' "
            + new_file_name
        )
        os.system("sed -i '27s@.*@  point 2.0 0.0 0.0" + "@' " + new_file_name)
        # Define optimization routine
        os.system("sed -i '32s@.*@  1 1 1 1 1 1 0" + "@' " + new_file_name)
        os.system("sed -i '33s@.*@  0 0 0 1 1 0 0" + "@' " + new_file_name)
        os.system("sed -i '34s@.*@  0 1 1" + "@' " + new_file_name)
        # Define executione commands
        os.system("sed -i '38s@.*@start_command@' " + new_file_name)
        os.system(
            "sed -i '40s@.*@readobs_point "
            + args["inbase"]
            + "/source_obs_"
            + str(ii)
            + ".dat"
            + "@' "
            + new_file_name
        )
        os.system(
            "sed -i '41s@.*@parprior "
            + args["inbase"]
            + "/prior_"
            + str(ii)
            + ".dat"
            + "@' "
            + new_file_name
        )
        SIEPOI_name = args["outbase"] + "/SIE_POI_" + str(ii)
        os.system("sed -i '43s@.*@optimize @' " + new_file_name)
        os.system("sed -i '44s@.*@resetopt_lens @' " + new_file_name)
        os.system("sed -i '45s@.*@calcein 2.0 @' " + new_file_name)
        os.system("sed -i '46s@.*@findimg @' " + new_file_name)
        
        # Not Working
        # subprocess.call(["./lensmodel", new_file_name])


if __name__ == "__main__":
    sl_sys_analysis()
