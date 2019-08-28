# — 1. modelled positions of the source (ys1, ys2)
# — 2. Einstein Radius of each lens (Re)
# — 3. modelled positions of the lensed images (x1_i, x2_i).
# — 4. modelled magnification of the lensed images (mu_i).
# — 5. modelled axis-ration and position angles of the lenses (ql and pa)

from __future__ import division
import os, sys, glob
import numpy as np
import pandas as pd
import xarrat as xr
import json


def make_map(ids, values, Nx, Ny):
    zip_array = sorted(zip(ids, values))
    ids = np.array([i for (i, j) in zip_array])
    values = np.array([j for (i, j) in zip_array])
    map_out = np.zeros([Nx, Ny])
    k = 0
    for j in range(Ny - 1, -1, -1):
        for i in range(Nx):
            map_out[j, i] = values[k]
            k += 1
    return map_out


args = {}
args["outdirstr"] = sys.argv[1]
args["los"] = sys.argv[2]
args["nimgs"] = sys.argv[3]
args["extkappa"] = sys.argv[4]
args["infile"] = sys.argv[5]

outdir = os.fsencode(args["outdirstr"])

resdir = []  # initialize output dictionary

# Run through files
files = profile_files = glob.glob(args["outdirstr"] + "fit6_*.chi")

init = 1
for indx, filename in enumerate(files):
    system_id = filename.split("_")[-1].split(".")[0]

    if init:
        df_deg = pd.read_csv(
            filename,
            sep=" ",
            skiprows=1,
            names=[
                "e",
                "shear",
                "chi2",
                "pos",
                "flux",
                "tdel",
                "gal",
                "plim",
                "crv",
                "ring",
            ],
            usecols=["e", "shear", "chi2"],
        )
        value_maps = np.zeros((df_deg["e"].values, df_deg["e"].values, len(files)))
        value_maps[:, :, indx] = make_map(
            df_deg.index.values,
            df_deg["chi2"].values,
            len(df_deg.index.values),
            len(df_deg.index.values),
        )
        init = 0
    else:
        df_deg = pd.read_csv(
            filename,
            sep=" ",
            skiprows=1,
            names=[
                "e",
                "shear",
                "chi2",
                "pos",
                "flux",
                "tdel",
                "gal",
                "plim",
                "crv",
                "ring",
            ],
            usecols=["e", "shear", "chi2"],
        )
        value_maps[:, :, indx] = make_map(
            df_deg.index.values,
            df_deg["chi2"].values,
            len(df_deg.index.values),
            len(df_deg.index.values),
        )


da = xr.DataArray(
    value_maps,
    coords=[df_deg["e"].values, df_deg["shear"].values],
    dims=["ellipticity", "shear"],
)
da.to_netcdf("./ellipticity_shear_deg.nc")

# Run through files
files = profile_files = glob.glob(args["outdirstr"] + "fitH0_*.chi")
init = 1
for filename in files:
    system_id = filename.split("_")[-1].split(".")[0]

    if init:
        df_H0 = pd.read_csv(
            filename,
            sep=" ",
            skiprows=1,
            names=[
                "h",
                "chi2_%s" % system_id,
                "pos",
                "flux",
                "tdel",
                "gal",
                "plim",
                "crv",
                "ring",
            ],
            index_col="h",
            usecols="chi2_%s" % system_id,
        )
        init = 0
    else:
        df_H0_new = pd.read_csv(
            filename,
            sep=" ",
            skiprows=1,
            names=[
                "h",
                "chi2_%s" % system_id,
                "pos",
                "flux",
                "tdel",
                "gal",
                "plim",
                "crv",
                "ring",
            ],
            index_col="h",
            usecols="chi2_%s" % system_id,
        )
        df_H0.join(df_H0_new)

df_H0.to_hdf("./EDA_H0_" + args["infile"], key="df_H0")
