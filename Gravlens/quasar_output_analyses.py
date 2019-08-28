# — 1. modelled positions of the source (ys1, ys2)
# — 2. Einstein Radius of each lens (Re)
# — 3. modelled positions of the lensed images (x1_i, x2_i).
# — 4. modelled magnification of the lensed images (mu_i).
# — 5. modelled axis-ration and position angles of the lenses (ql and pa)

from __future__ import division
import os, sys, glob
import numpy as np
import pandas as pd
import xarray as xr
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
            skiprows=3,
            names=[
                "e",
                "shear",
                "chi2",
            ],
            usecols=["e", "shear", "chi2"],
        )
        cell_nr = len(df_deg.index.values)
        if (cell_nr == 50*50) and \
           (len(np.unique(df_deg["e"].values)) == 50) and \
           (len(np.unique(df_deg["shear"].values)) == 50):
            print(filename)
            grid_cells_x = int(np.sqrt(len(df_deg.index.values)))
            grid_cells_y = int(np.sqrt(len(df_deg.index.values)))
            ellipticity_coord = np.unique(df_deg["e"].values)
            shear_coord = np.unique(df_deg["shear"].values)
            value_maps = np.zeros((grid_cells_x, grid_cells_y, len(files)))
            value_maps[:, :, indx] = make_map(
                df_deg.index.values,
                df_deg["chi2"].values,
                grid_cells_x,
                grid_cells_y,
            )
            init = 0
        else:
            continue
    else:
        df_deg = pd.read_csv(
            filename,
            sep=" ",
            skiprows=3,
            names=[
                "e",
                "shear",
                "chi2",
            ],
            usecols=["e", "shear", "chi2"],
        )
        cell_nr = len(df_deg.index.values)
        if (cell_nr == 50*50) and \
           (len(np.unique(df_deg["e"].values)) == 50) and \
           (len(np.unique(df_deg["shear"].values)) == 50):
            print(filename)
            value_maps[:, :, indx] = make_map(
                df_deg.index.values,
                df_deg["chi2"].values,
                grid_cells_x,
                grid_cells_y,
            )
        else:
            continue

da = xr.DataArray(
    value_maps,
    coords=[ellipticity_coord, shear_coord, np.arange(len(files))],
    dims=["ellipticity", "shear", "lensID"],
)
da.to_netcdf("./check_ellipticity_shear_deg.nc")

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
            ],
            usecols=["h", "chi2_%s" % system_id],
            index_col="h",
        )
        if len(df_H0.index.values) == 101:
            print(filename)
            init = 0
    else:
        df_H0_new = pd.read_csv(
            filename,
            sep=" ",
            skiprows=1,
            names=[
                "chi2_%s" % system_id,
            ],
            usecols=["chi2_%s" % system_id],
        )
        if len(df_H0_new.index.values) == 101:
            print(filename)
            df_H0.join(df_H0_new)

df_H0.to_hdf("./check_H0_" + args["infile"], key="df_H0")
