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

args = {}
args["outdirstr"] = sys.argv[1]
args["los"] = sys.argv[2]
args["nimgs"] = sys.argv[3]
args["extkappa"] = sys.argv[4]
args["infile"] = sys.argv[5]

outdir = os.fsencode(args["outdirstr"])

resdir = []  # initialize output dictionary

# Run through files
files = profile_files = glob.glob(args["outdirstr"] + "fit6_*.best")

for filename in files:
    system_id = filename.split("_")[-1].split(".")[0]
