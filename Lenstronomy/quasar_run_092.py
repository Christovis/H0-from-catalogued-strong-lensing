import sys
import numpy as np
import corner
import lenstronomy
from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.Sampling.parameters import Param
import lenstronomy.Util.param_util as param_util
from lenstronomy.Analysis.lens_analysis import LensAnalysis
from lenstronomy.Cosmo.lens_cosmo import LensCosmo

# import fastell4py
import time
from glob import glob
import json
import h5py
from mpi4py import MPI
from mpi_errchk import mpi_errchk

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
        args["version"] = sys.argv[4]
        args["dt_sigma"] = float(sys.argv[5])
        args["image_amps_sigma"] = float(sys.argv[6])
        args["flux_ratio_errors"] = float(sys.argv[7])
        args["astrometry_sigma"] = float(sys.argv[8])

    args = comm.bcast(args)
    # Organize devision of strong lensing systems
    with open(args["infile"], "r") as myfile:
        limg_data = myfile.read()
    systems = json.loads(limg_data)
    sys_nr_per_proc = int(len(systems) / comm_size)
    print("comm_rank", comm_rank)
    start_sys = sys_nr_per_proc * comm_rank
    end_sys = sys_nr_per_proc * (comm_rank + 1)
    print(start_sys, end_sys)
    with open("../lens_catalogs_sie_only.json", "r") as myfile:
        limg_data = myfile.read()
    systems_prior = json.loads(limg_data)

    if comm_rank == 0:
        print("Each process will have %d systems" % sys_nr_per_proc)
        print("That should take app. %f min." % (sys_nr_per_proc * 20))

    results = {"gamma": [], "phi_ext": [], "gamma_ext": [], "theta_E": [], "D_dt": []}
    for ii in range(len(systems))[(start_sys + 2) : end_sys]:
        system = systems[ii]
        system_prior = systems_prior[ii]
        print("Analysing system ID: %d" % ii)
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
            ximg[jj] = system["ximg"][jj]  # [arcsec]
            yimg[jj] = system["yimg"][jj]  # [arcsec]
            delay[jj] = system["delay"][jj]  # [days]
            image_amps[jj] = system["mags"][jj]  # [linear units or magnitudes]
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
        # keyword list with all the data elements.
        kwargs_data_joint = {
            "time_delays_measured": d_dt,
            "time_delays_uncertainties": d_dt_sigma,
            "flux_ratios": flux_ratios,
            "flux_ratio_errors": flux_ratio_errors,
        }

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
        # error
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
        # lower limit
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
        # upper limit
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
        kwargs_lens_init.append({"gamma_ext": 0.05, "psi_ext": 0.0})
        kwargs_lens_sigma.append({"gamma_ext": 0.05, "psi_ext": np.pi})
        kwargs_lower_lens.append({"gamma_ext": 0, "psi_ext": -np.pi})
        kwargs_upper_lens.append({"gamma_ext": 0.3, "psi_ext": np.pi})
        # combined lens model
        lens_params = [
            kwargs_lens_init,
            kwargs_lens_sigma,
            fixed_lens,
            kwargs_lower_lens,
            kwargs_upper_lens,
        ]

        # image position parameters
        point_source_list = ["LENSED_POSITION"]
        # we fix the image position coordinates
        fixed_ps = [{}]
        # the initial guess for the appearing image positions is:
        # at the image position.
        kwargs_ps_init = [{"ra_image": ximg, "dec_image": yimg}]
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
            {"ra_image": 10 * np.ones(len(ximg)), "dec_image": 10 * np.ones(len(ximg))}
        ]

        ps_params = [
            kwargs_ps_init,
            kwargs_ps_sigma,
            fixed_ps,
            kwargs_lower_ps,
            kwargs_upper_ps,
        ]

        # quasar source size
        fixed_special["source_size"] = source_size_arcsec
        kwargs_special_init["source_size"] = source_size_arcsec
        kwargs_special_sigma["source_size"] = source_size_arcsec
        kwargs_lower_special["source_size"] = 0.0001

        # Time-delay distance
        kwargs_special_init["D_dt"] = 5000
        kwargs_special_sigma["D_dt"] = 3000
        kwargs_lower_special["D_dt"] = 0
        kwargs_upper_special["D_dt"] = 10000

        special_params = [
            kwargs_special_init,
            kwargs_special_sigma,
            fixed_special,
            kwargs_lower_special,
            kwargs_upper_special,
        ]

        # combined parameter settings
        kwargs_params = {
            "lens_model": lens_params,
            "point_source_model": ps_params,
            "special": special_params,
        }

        # our model choices
        kwargs_model = {
            "lens_model_list": lens_model_list,
            "point_source_model_list": point_source_list,
        }

        # setup options for likelihood and parameter sampling
        kwargs_constraints = {
            "num_point_source_list": [int(args["nimgs"])],
            # any proposed lens model must satisfy the image positions
            # appearing at the position of the point sources being sampeld
            "solver_type": "PROFILE_SHEAR",
            "Ddt_sampling": time_delay_likelihood,
            # sampling of the time-delay distance
            # explicit modelling of the astrometric imperfection of
            # the point source positions
            "point_source_offset": True,
        }

        # explicit sampling of finite source size parameter
        # (only use when source_type='GAUSSIAN' or 'TORUS')
        if (
            kwargs_flux_compute["source_type"] in ["GAUSSIAN", "TORUS"]
            and flux_ratio_likelihood is True
        ):
            kwargs_constraints["source_size"] = True

        # e.g. power-law mass slope of the main deflector
        # [[index_model, 'param_name', mean, 1-sigma error], [...], ...]
        prior_lens = [[0, "gamma", 2, 0.1]]
        prior_special = []

        kwargs_likelihood = {
            "position_uncertainty": args["astrometry_sigma"],
            "source_position_likelihood": True,
            "image_position_likelihood": True,
            "time_delay_likelihood": True,
            "flux_ratio_likelihood": True,
            "kwargs_flux_compute": kwargs_flux_compute,
            "prior_lens": prior_lens,
            "prior_special": prior_special,
            "check_solver": True,
            "solver_tolerance": 0.001,
            "check_bounds": True,
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
            ["PSO", {"sigma_scale": 1.0, "n_particles": 100, "n_iterations": 100}]
        ]

        chain_list_pso = fitting_seq.fit_sequence(fitting_kwargs_list)
        kwargs_result = fitting_seq.best_fit()
        kwargs_result = fitting_seq.best_fit(bijective=True)
        args_result = fitting_seq.param_class.kwargs2args(**kwargs_result)
        logL, _ = fitting_seq.likelihoodModule.logL(args_result, verbose=True)

        # and now we run the MCMC
        fitting_kwargs_list = [
            [
                "MCMC",
                {"n_burn": 400, "n_run": 600, "walkerRatio": 10, "sigma_scale": 0.1},
            ]
        ]
        chain_list_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
        kwargs_result = fitting_seq.best_fit()
        # print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
        # print("parameters in order: ", param_mcmc)
        print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])

        param = Param(
            kwargs_model,
            fixed_lens,
            kwargs_fixed_ps=fixed_ps,
            kwargs_fixed_special=fixed_special,
            kwargs_lens_init=kwargs_result["kwargs_lens"],
            **kwargs_constraints,
        )
        # the number of non-linear parameters and their names #
        num_param, param_list = param.num_param()

        lensModel = LensModel(kwargs_model["lens_model_list"])
        lensModelExtensions = LensModelExtensions(lensModel=lensModel)

        mcmc_new_list = []
        labels_new = [
            r"$\theta_E$",
            r"$\gamma$",
            r"$\phi_{lens}$",
            r"$q$",
            r"$\phi_{ext}$",
            r"$\gamma_{ext}$",
        ]
        if flux_ratio_likelihood is True:
            labels_new.extend(["B/A", "C/A", "D/A"])
        if kwargs_constraints.get("source_size", False) is True:
            if "source_size" not in fixed_special:
                labels_new.append("source size")
        if time_delay_likelihood is True:
            labels_new.append(r"$D_{dt}$")

        for i in range(len(samples_mcmc)):
            kwargs_out = param.args2kwargs(samples_mcmc[i])
            kwargs_lens_out, kwargs_special_out, kwargs_ps_out = (
                kwargs_out["kwargs_lens"],
                kwargs_out["kwargs_special"],
                kwargs_out["kwargs_ps"],
            )

            # compute 'real' image position adding potential astrometric shifts
            x_pos = kwargs_ps_out[0]["ra_image"]
            y_pos = kwargs_ps_out[0]["dec_image"]

            # extract quantities of the main deflector
            theta_E = kwargs_lens_out[0]["theta_E"]
            gamma = kwargs_lens_out[0]["gamma"]
            e1, e2 = kwargs_lens_out[0]["e1"], kwargs_lens_out[0]["e2"]
            phi, q = param_util.ellipticity2phi_q(e1, e2)
            phi_ext, gamma_ext = (
                kwargs_lens_out[1]["psi_ext"] % np.pi,
                kwargs_lens_out[1]["gamma_ext"],
            )
            new_chain = [theta_E, gamma, phi, q, phi_ext, gamma_ext]
            if flux_ratio_likelihood is True:
                mag = lensModel.magnification(x_pos, y_pos, kwargs_lens_out)
                flux_ratio_fit = mag[1:] / mag[0]
                new_chain.extend(
                    [flux_ratio_fit[0], flux_ratio_fit[1], flux_ratio_fit[2]]
                )
            if (
                kwargs_constraints.get("source_size", False) is True
                and "source_size" not in fixed_special
            ):
                source_size = kwargs_special_out["source_size"]
                new_chain.append(source_size)
            if time_delay_likelihood is True:
                D_dt = kwargs_special_out["D_dt"]
                new_chain.append(D_dt)

            # source_size = np.random.uniform(high=1, low=0)
            mcmc_new_list.append(np.array(new_chain))

        # plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
        # and here the predicted angular diameter distance from a
        # default cosmology (attention for experimenter bias!)
        cosmo = FlatLambdaCDM(H0=71, Om0=0.3089, Ob0=0.0)
        lensCosmo = LensCosmo(cosmo=cosmo, z_lens=z_lens, z_source=z_source)
        gamma = np.mean(gamma)
        phi_ext = np.mean(phi_ext)
        gamma_ext = np.mean(gamma_ext)
        theta_E = np.mean(theta_E)
        D_dt = np.mean(D_dt)
        results["gamma"].append(gamma)
        results["phi_ext"].append(phi_ext)
        results["gamma_ext"].append(gamma_ext)
        results["theta_E"].append(theta_E)
        results["D_dt"].append(lensCosmo.D_dt)

    with open(
        "./quasars_%s_nimgs_%s_%s.json" % (args["los"], args["nimgs"], args["version"]),
        "w",
    ) as fout:
        json.dump(results, fout)


if __name__ == "__main__":
    sl_sys_analysis()
