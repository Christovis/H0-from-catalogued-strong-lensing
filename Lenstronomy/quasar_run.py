import sys
import numpy as np
import corner
import lenstronomy
from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.Sampling.parameters import Param
from lenstronomy.Analysis.lens_analysis import LensAnalysis
from lenstronomy.Cosmo.lens_cosmo import LensCosmo

# import fastell4py
import time
from glob import glob
import json
from mpi4py import MPI
from mpi_errchk import mpi_errchk
import h5py

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
        args["nimgs"] = sys.argv[2]
        args["los"] = sys.argv[3]
        args["dt_sigma"] = float(sys.argv[4])
        args["image_amps_sigma"] = float(sys.argv[5])
        args["flux_ratio_errors"] = float(sys.argv[6])
        args["astrometry_sigma"] = float(sys.argv[7])

    args = comm.bcast(args)
    # Organize devision of strong lensing systems
    with open(args["infile"], "r") as myfile:
        limg_data = myfile.read()
    systems = json.loads(limg_data)
    sys_nr_per_proc = int(len(systems) / comm_size)
    print('comm_rank', comm_rank)
    start_sys = sys_nr_per_proc * comm_rank
    end_sys = sys_nr_per_proc * (comm_rank + 1)
    print(start_sys, end_sys)
    with open("../lens_catalogs_sie_only.json", "r") as myfile:
        limg_data = myfile.read()
    systems_prior = json.loads(limg_data)

    if comm_rank == 0:
        print("Each process will have %d systems" % sys_nr_per_proc)
        print("That should take app. %f min." % (sys_nr_per_proc * 20))

    results = {"D_dt": []}
    for ii in range(len(systems))[(start_sys + 2) : end_sys]:
        system = systems[ii]
        system_prior = systems_prior[ii]
        print("Analysing system ID: %d" % system["losID"])
        print(system)
        # the data set is
        z_lens = system_prior["zl"]
        z_source = 2.0

        # multiple images properties
        ximg = np.zeros(system["nimgs"])
        yimg = np.zeros(system["nimgs"])
        delay = np.zeros(system["nimgs"])
        image_amps = np.zeros(system["nimgs"])
        for jj in range(system["nimgs"]):
            ximg[jj] = system["ximg"][jj]        #[arcsec]
            yimg[jj] = system["yimg"][jj]        #[arcsec]
            delay[jj] = system["delay"][jj]      #[days]
            image_amps[jj] = system["mags"][jj]  #[linear units or magnitudes]
        # sort by arrival time
        index_sort = np.argsort(delay)
        ximg = ximg[index_sort]
        yimg = yimg[index_sort]
        delay = delay[index_sort]
        image_amps = image_amps[index_sort]
        d_dt = delay[1:] - delay[0]
        # measurement uncertainties
        d_dt_sigma = np.ones(system["nimgs"] - 1) * args["dt_sigma"]
        image_amps_sigma = np.ones(system["nimgs"]) * args["image_amps_sigma"]
        flux_ratios = image_amps[1:] - image_amps[0]
        flux_ratio_errors = np.ones(system["nimgs"] - 1) * args["flux_ratio_errors"]

        # lens model choices
        lens_model_list = ["SPEP", "SHEAR"]
        # first choice: SPEP
        fixed_lens = []
        kwargs_lens_init = []
        kwargs_lens_sigma = []
        kwargs_lower_lens = []
        kwargs_upper_lens = []
        fixed_lens.append({})
        kwargs_lens_init.append(
            {
                "theta_E": 1.0,
                "gamma": 2,
                "center_x": 0,
                "center_y": 0,
                "e1": 0,
                "e2": 0.0,
            }
        )
        kwargs_lens_sigma.append(
            {
                "theta_E": 0.2,
                "e1": 0.1,
                "e2": 0.1,
                "gamma": 0.1,
                "center_x": 0.1,
                "center_y": 0.1,
            }
        )
        kwargs_lower_lens.append(
            {
                "theta_E": 0.01,
                "e1": -0.5,
                "e2": -0.5,
                "gamma": 1.5,
                "center_x": -10,
                "center_y": -10,
            }
        )
        kwargs_upper_lens.append(
            {
                "theta_E": 10,
                "e1": 0.5,
                "e2": 0.5,
                "gamma": 2.5,
                "center_x": 10,
                "center_y": 10,
            }
        )
        # second choice: SHEAR
        fixed_lens.append({"ra_0": 0, "dec_0": 0})
        kwargs_lens_init.append({"e1": 0.0, "e2": 0.0})
        kwargs_lens_sigma.append({"e1": 0.1, "e2": 0.1})
        kwargs_lower_lens.append({"e1": -0.2, "e2": -0.2})
        kwargs_upper_lens.append({"e1": 0.2, "e2": 0.2})
        lens_params = [
            kwargs_lens_init,
            kwargs_lens_sigma,
            fixed_lens,
            kwargs_lower_lens,
            kwargs_upper_lens,
        ]

        point_source_list = ["LENSED_POSITION"]
        fixed_ps = [
            {"ra_image": ximg, "dec_image": yimg}
        ]
        kwargs_ps_init = fixed_ps
        # let some freedome in how well the actual image positions are
        # matching those given by the data (indicated as 'ra_image', 'dec_image'
        # and held fixed while fitting)
        kwargs_ps_sigma = [
            {
                "ra_image": 0.01 * np.ones(len(ximg)),
                "dec_image": 0.01 * np.ones(len(ximg)),
            }
        ]
        kwargs_lower_ps = [
            {
                "ra_image": -10 * np.ones(len(ximg)),
                "dec_image": -10 * np.ones(len(ximg)),
            }
        ]
        kwargs_upper_ps = [
            {"ra_image": 10 * np.ones(len(ximg)),
             "dec_image": 10 * np.ones(len(ximg))}
        ]

        ps_params = [
            kwargs_ps_init,
            kwargs_ps_sigma,
            fixed_ps,
            kwargs_lower_ps,
            kwargs_upper_ps,
        ]

        fixed_cosmo = {}
        kwargs_cosmo_init = {
            "D_dt": 5000,
            "delta_x_image": np.zeros_like(ximg),
            "delta_y_image": np.zeros_like(ximg),
        }
        kwargs_cosmo_sigma = {
            "D_dt": 10000,
            "delta_x_image": np.ones_like(ximg) * args["astrometry_sigma"],
            "delta_y_image": np.ones_like(ximg) * args["astrometry_sigma"],
        }
        kwargs_lower_cosmo = {
            "D_dt": 0,
            "delta_x_image": np.ones_like(ximg) * (-1),
            "delta_y_image": np.ones_like(ximg) * (-1),
        }
        kwargs_upper_cosmo = {
            "D_dt": 10000,
            "delta_x_image": np.ones_like(ximg) * (1),
            "delta_y_image": np.ones_like(ximg) * (1),
        }
        cosmo_params = [
            kwargs_cosmo_init,
            kwargs_cosmo_sigma,
            fixed_cosmo,
            kwargs_lower_cosmo,
            kwargs_upper_cosmo,
        ]

        kwargs_params = {
            "lens_model": lens_params,
            "point_source_model": ps_params,
            "cosmography": cosmo_params,
        }

        # setup options for likelihood and parameter sampling
        kwargs_constraints = {
            "num_point_source_list": [int(args["nimgs"])],
            # any proposed lens model must satisfy the image positions
            # appearing at the position of the point sources being sampeld
            "solver_type": "PROFILE_SHEAR",
            "cosmo_type": "D_dt",
            # sampling of the time-delay distance
            # explicit modelling of the astrometric imperfection of
            # the point source positions
            "point_source_offset": True,
        }
        kwargs_likelihood = {
            "check_bounds": True,
            "point_source_likelihood": True,
            "position_uncertainty": args["astrometry_sigma"],
            "check_solver": True,
            "solver_tolerance": 0.001,
            "time_delay_likelihood": True,
            "image_likelihood": False,
            # this needs to be explicitly given when not having imaging data
            "flux_ratio_likelihood": True,  # enables the flux ratio likelihood
        }
        kwargs_data_joint = {
            "time_delays_measured": d_dt,
            "time_delays_uncertainties": d_dt_sigma,
            "flux_ratios": flux_ratios,
            "flux_ratio_errors": flux_ratio_errors,
        }
        kwargs_model = {
            "lens_model_list": lens_model_list,
            "point_source_model_list": point_source_list,
        }

        mpi = False  # MPI possible, but not supported through that notebook.
        fitting_seq = FittingSequence(
            kwargs_data_joint,
            kwargs_model,
            kwargs_constraints,
            kwargs_likelihood,
            kwargs_params,
        )
        fitting_kwargs_list = [
            # ['update_settings', {'lens_add_fixed': [[0, ['gamma']]]}],
            ["PSO", {"sigma_scale": 1.0, "n_particles": 100, "n_iterations": 100}]
        ]

        chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(
            fitting_kwargs_list
        )
        lens_result, source_result, lens_light_result, ps_result, cosmo_result = (
            fitting_seq.best_fit()
        )

        # and now we run the MCMC
        fitting_kwargs_list = [
            ["PSO", {"sigma_scale": 0.1, "n_particles": 100, "n_iterations": 100}],
            [
                "MCMC",
                {"n_burn": 200, "n_run": 200, "walkerRatio": 10, "sigma_scale": 0.1},
            ],
        ]
        chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(
            fitting_kwargs_list
        )
        lens_result, source_result, lens_light_result, ps_result, cosmo_result = (
            fitting_seq.best_fit()
        )
        # print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
        # print("parameters in order: ", param_mcmc)
        print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])

        param = Param(
            kwargs_model,
            fixed_lens,
            kwargs_fixed_ps=fixed_ps,
            kwargs_fixed_cosmo=fixed_cosmo,
            kwargs_lens_init=lens_result,
            **kwargs_constraints,
        )
        # the number of non-linear parameters and their names #
        num_param, param_list = param.num_param()

        lensAnalysis = LensAnalysis(kwargs_model)

        mcmc_new_list = []
        labels_new = [
            r"$\gamma$",
            r"$\phi_{ext}$",
            r"$\gamma_{ext}$",
            r"$D_{\Delta t}$",
        ]
        for i in range(len(samples_mcmc)):
            # transform the parameter position of the MCMC chain in a
            # lenstronomy convention with keyword arguments #
            kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out, kwargs_cosmo = param.args2kwargs(
                samples_mcmc[i]
            )
            D_dt = kwargs_cosmo["D_dt"]
            gamma = kwargs_lens_out[0]["gamma"]
            phi_ext, gamma_ext = lensAnalysis._lensModelExtensions.external_shear(
                kwargs_lens_out
            )
            mcmc_new_list.append([gamma, phi_ext, gamma_ext, D_dt])

        # plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
        # and here the predicted angular diameter distance from a
        # default cosmology (attention for experimenter bias!)
        lensCosmo = LensCosmo(z_lens=z_lens, z_source=z_source)
        results["D_dt"].append(D_dt)

    hf = h5py.File(
        "H0_"
        + args["los"]
        + "_nimgs"
        + str(args["nimgs"])
        + "_cpu"
        + str(comm_rank)
        + ".hdf5",
        "w",
    )
    hf.create_dataset("D_dt", data=results["D_dt"])
    hf.close()


if __name__ == "__main__":
    sl_sys_analysis()
