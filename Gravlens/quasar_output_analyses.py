# — 1. modelled positions of the source (ys1, ys2)
# — 2. Einstein Radius of each lens (Re)
# — 3. modelled positions of the lensed images (x1_i, x2_i).
# — 4. modelled magnification of the lensed images (mu_i).
# — 5. modelled axis-ration and position angles of the lenses (ql and pa)

from __future__ import division
import os, sys, glob
import numpy as np
import pandas as pd
import json


def _make_map(ids, values, Nx, Ny):
    zip_array = sorted(zip(ids, values))
    ids = np.array([i for (i,j) in zip_array])
    values = np.array([j for (i,j) in zip_array])
    map_out = np.zeros([Nx, Ny])
    k = 0
    for j in range(Ny-1, -1, -1):
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

value_maps = []
for filename in files:
    system_id = filename.split("_")[-1].split(".")[0]
    df_deg = pd.read_csv(
        filename, sep=" ", skiprows=1,
        names=['e', 'shear', 'chi2', 'pos', 'flux', 'tdel', 'gal', 'plim', 'crv', 'ring'],
        usecols=['e', 'shear', 'chi2'],
    )
    grid_cells = len(df_deg.index.values)
    value_map = _make_map(
        grid_cells,
        df_deg["chi2"].values,
        df_deg["e"].values,
        df_deg["shear"].values,
    )
    value_maps.append(value_map)

with open("./EDA_ellip_shear_" + args["infile"], 'w') as f:
    json.dump(value_maps, f)

# Run through files
files = profile_files = glob.glob(args["outdirstr"] + "fitH0_*.chi")
init = 0
for filename in files:
    system_id = filename.split("_")[-1].split(".")[0]

    if init:
        df_H0 = pd.read_csv(
            filename, sep=" ", skiprows=1,
            names=['h', 'chi2', 'pos', 'flux', 'tdel', 'gal', 'plim', 'crv', 'ring'],
            usecols=['h', 'chi2_%s' % system_id]
        )
    else:
        df_H0 = pd.read_csv(
            filename, sep=" ", skiprows=1,
            names=['h', 'chi2', 'pos', 'flux', 'tdel', 'gal', 'plim', 'crv', 'ring'],
            usecols=['chi2_%s' % system_id]
        )

df_H0.to_hdf("./EDA_H0_" + args["infile"], key='df_H0')
