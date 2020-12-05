import sys
import numpy as np
import corner

from astropy.cosmology import FlatLambdaCDM

import lenstronomy
from lenstronomy.Util import constants
import lenstronomy.Util.param_util as param_util
from lenstronomy.Sampling.parameters import Param
from lenstronomy.Cosmo.lens_cosmo import LensCosmo
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.Analysis.lens_analysis import LensAnalysis
from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.LensModel.lens_model_extensions import LensModelExtensions
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver

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

    source_size_pc = 10.0
    window_size = 0.1  # units of arcseconds
    grid_number = 100  # supersampled window (per axis)
    z_source = 2.0
    cosmo = FlatLambdaCDM(H0=71, Om0=0.3089, Ob0=0.0)

    results = {"gamma": [], "phi_ext": [], "gamma_ext": [], "theta_E": [], "D_dt": []}
    for ii in range(len(systems))[(start_sys + 2) : end_sys]:
        system = systems[ii]
        system_prior = systems_prior[ii]
        print("Analysing system ID: %d" % ii)

        # the data set is
        z_lens = system_prior["zl"]
        lensCosmo = LensCosmo(cosmo=cosmo, z_lens=z_lens, z_source=z_source)
        # convert units of pc into arcseconds
        D_s = lensCosmo.D_s
        source_size_arcsec = source_size_pc / 10 ** 6 / D_s / constants.arcsec
        print("The source size in arcsec init = %.4f" % source_size_arcsec) #0.0012

        # multiple images properties
        ximg = np.zeros(system["nimgs"])
        yimg = np.zeros(system["nimgs"])
        t_days = np.zeros(system["nimgs"])
        image_amps = np.zeros(system["nimgs"])
        for jj in range(system["nimgs"]):
            ximg[jj] = system["ximg"][jj]  # [arcsec]
            yimg[jj] = system["yimg"][jj]  # [arcsec]
            t_days[jj] = system["delay"][jj]  # [days]
            image_amps[jj] = system["mags"][jj]  # [linear units or magnitudes]
        # sort by arrival time
        index_sort = np.argsort(t_days)
        ximg = ximg[index_sort]  # relative RA (arc seconds)
        yimg = yimg[index_sort]  # relative DEC (arc seconds)
        image_amps = np.abs(image_amps[index_sort])
        t_days = t_days[index_sort]
        d_dt = t_days[1:] - t_days[0]

        # measurement uncertainties
        astrometry_sigma = args["astrometry_sigma"]
        ximg_measured = ximg + np.random.normal(0, astrometry_sigma, system["nimgs"])
        yimg_measured = yimg + np.random.normal(0, astrometry_sigma, system["nimgs"])
        image_amps_sigma = np.ones(system["nimgs"]) * args["image_amps_sigma"]
        flux_ratios = image_amps[1:] - image_amps[0]
        flux_ratio_errors = np.ones(system["nimgs"] - 1) * args["flux_ratio_errors"]
        flux_ratios_measured = flux_ratios + np.random.normal(0, flux_ratio_errors)
        d_dt_sigma = np.ones(system["nimgs"] - 1) * args["dt_sigma"]
        d_dt_measured = d_dt + np.random.normal(0, d_dt_sigma)

        kwargs_data_joint = {
            "time_delays_measured": d_dt_measured,
            "time_delays_uncertainties": d_dt_sigma,
            "flux_ratios": flux_ratios_measured,
            "flux_ratio_errors": flux_ratio_errors,
            "ra_image_list": [ximg_measured],
            "dec_image_list": [yimg_measured],
        }

        # lens model choices
        lens_model_list = ["SPEMD", "SHEAR_GAMMA_PSI"]

        # 1. layer: primary SPEP
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
        # 2nd layer: external SHEAR
        fixed_lens.append({"ra_0": 0, "dec_0": 0})
        kwargs_lens_init.append({"gamma_ext": 0.05, "psi_ext": 0.0})
        kwargs_lens_sigma.append({"gamma_ext": 0.05, "psi_ext": np.pi})
        kwargs_lower_lens.append({"gamma_ext": 0, "psi_ext": -np.pi})
        kwargs_upper_lens.append({"gamma_ext": 0.3, "psi_ext": np.pi})
        
        # 3rd layer: external CONVERGENCE
		kwargs_lens_init.append({'kappa_ext': 0.12})
		kwargs_lens_sigma.append({'kappa_ext': 0.06})
		kwargs_lower_lens.append({'kappa_ext': 0.0})
		kwargs_upper_lens.append({'kappa_ext': 0.3})
        
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
        fixed_special = {}
        kwargs_special_init = {}
        kwargs_special_sigma = {}
        kwargs_lower_special = {}
        kwargs_upper_special = {}

        fixed_special["source_size"] = source_size_arcsec
        kwargs_special_init["source_size"] = source_size_arcsec
        kwargs_special_sigma["source_size"] = source_size_arcsec
        kwargs_lower_special["source_size"] = 0.0001
        kwargs_upper_special["source_size"] = 1

        # Time-delay distance
        kwargs_special_init["D_dt"] = 4300   # corresponds to H0 ~ 70
        kwargs_special_sigma["D_dt"] = 3000
        kwargs_lower_special["D_dt"] = 2500  # corresponds to H0 ~ 120
        kwargs_upper_special["D_dt"] = 14000 # corresponds to H0 ~ 20

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
        lensModel = LensModel(kwargs_model["lens_model_list"])
        lensModelExtensions = LensModelExtensions(lensModel=lensModel)
        lensEquationSolver = LensEquationSolver(lensModel=lensModel)

        # setup options for likelihood and parameter sampling
        time_delay_likelihood = True
        flux_ratio_likelihood = True
        image_position_likelihood = True
        kwargs_flux_compute = {
            "source_type": "INF",
            "window_size": window_size,
            "grid_number": grid_number,
        }

        kwargs_constraints = {
            "num_point_source_list": [int(args["nimgs"])],
            # any proposed lens model must satisfy the image positions
            # appearing at the position of the point sources being sampeld
            # "solver_type": "PROFILE_SHEAR",
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

        fitting_seq = FittingSequence(
            kwargs_data_joint,
            kwargs_model,
            kwargs_constraints,
            kwargs_likelihood,
            kwargs_params,
        )
        fitting_kwargs_list = [
            ["PSO", {"sigma_scale": 1.0, "n_particles": 200, "n_iterations": 500}]
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
            if flux_ratio_likelihood is True:
                mag = lensModel.magnification(x_pos, y_pos, kwargs_lens_out)
                flux_ratio_fit = mag[1:] / mag[0]
            if (
                kwargs_constraints.get("source_size", False) is True
                and "source_size" not in fixed_special
            ):
                source_size = kwargs_special_out["source_size"]
            if time_delay_likelihood is True:
                D_dt = kwargs_special_out["D_dt"]

        # and here the predicted angular diameter distance from a
        # default cosmology (attention for experimenter bias!)
        gamma = np.median(gamma)
        phi_ext = np.median(phi_ext)
        gamma_ext = np.median(gamma_ext)
        theta_E = np.median(theta_E)
        D_dt = np.median(D_dt)
        results["gamma"].append(gamma)
        results["phi_ext"].append(phi_ext)
        results["gamma_ext"].append(gamma_ext)
        results["theta_E"].append(theta_E)
        results["H0"].append(c_light / D_dt)

    with open("./results_%s" % (args["infile"]), "w") as fout:
        json.dump(results, fout)


if __name__ == "__main__":
    sl_sys_analysis()
