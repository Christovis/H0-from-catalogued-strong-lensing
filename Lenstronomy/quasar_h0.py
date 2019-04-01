import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lensutils as lu
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.PointSource.point_source import PointSource


# lensing quantities
kwargs_shear = {
    'e1': -0.04,
    'e2': -0.01}  # shear values to the source plane
kwargs_spemd = {
    'theta_E': 1.66,
    'gamma': 1.98,
    'center_x': 0.0,
    'center_y': 0.0,
    'e1': 0.05,
    'e2': 0.05}  # parameters of the deflector lens model

# the lens model is a supperposition of an elliptical lens model with external shear
lens_model_list = ['SPEP', 'SHEAR']
kwargs_lens = [kwargs_spemd, kwargs_shear]
lens_model_class = LensModel(lens_model_list=lens_model_list)

# choice of source type
source_type = 'POINT'
source_x = 0.0
source_y = 0.0

indir = $1
with open(args["infile"], 'r') as myfile:
    limg_data = myfile.read()
systems = json.loads(limg_data)
# lesnes = list of all 1000 lenses by Nan
lenses = []
for ii in range(len(systems)):
    system = systems[ii]
    (system["ximg"][jj], system["yimg"][jj],
            system["mags"][jj]
            system["delay"][jj],
    lenses.append(
        lu.LenstronomyLens(
            name=nimgs+los+ii, longname=nimgs+los+ii,
            zlens=0.50175, zsource=2.0,
        )

# lensed image positions (solution of the lens equation) #
point_source_model_list = ['LENSED_POSITION']
pointSource = PointSource(
    point_source_type_list=point_source_model_list,
    lensModel=lensModel,
    fixed_magnification_list=[False])
kwargs_ps = [{
    'ra_image': theta_ra,
    'dec_image': theta_dec,
    'point_amp': np.abs(mag)*30 # the magnification of the point source images
    }]
# return image positions and amplitudes #
x_pos, y_pos = pointSource.image_position(
    kwargs_ps=kwargs_ps,
    kwargs_lens=kwargs_lens)
point_amp = pointSource.image_amplitude(
    kwargs_ps=kwargs_ps, kwargs_lens=kwargs_lens)


# MCMC parameters --------------------------------------------------------------------
compute_individual_samples = True  # if true the cosmo. posteriors are sampled for each lens
saveresults = True  # if true MCMC chains computed with emcee are saved
outdir = "samples"  # output directory
nwalkers = 6  #32
nsamples = 7000  #20000
cosmology = "FLCDM"

# MCMC sampling -----------------------------------------------------------------
display("Sampling cosmological parameters")
savedir = os.path.join(outdir, cosmology, "%ix%i" % (nwalkers, nsamples))
if not os.path.exists(savedir):
    display("Creating directory %s" % savedir)
    os.makedirs(savedir)    

combined_samples = lu.sample_params(
        lenses, cosmology=cosmology, nwalkers=nwalkers,
        nsamples=nsamples, save=saveresults, 
        filepath=os.path.join(savedir, "sample_%s.pkl" 
                                        % ("+".join([l.name for l in lenses])))
        )


