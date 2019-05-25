# — 1. modelled positions of the source (ys1, ys2)
# — 2. Einstein Radius of each lens (Re)
# — 3. modelled positions of the lensed images (x1_i, x2_i).
# — 4. modelled magnification of the lensed images (mu_i).
# — 5. modelled axis-ration and position angles of the lenses (ql and pa) 

from __future__ import division
import os, sys
import numpy as np
import pandas as pd
import json

args = {}
args["outdirstr"] = sys.argv[1]
args["los"] = sys.argv[2]
args["nimgs"] = sys.argv[3]
outdir = os.fsencode(args["outdirstr"])

mass_scale = []
lens_pos = []
ellipticity = []
position_angle = []
shear = []
shear_angle = []
hubble_const = []
resdir = []


# Run through files
for file in os.listdir(outdir):
    filename = os.fsdecode(file)
    system_id = filename.split('_')[1]

    with open(args["outdirstr"]+filename, "r") as f:
        lines = f.readlines()

        # Run through lines in file
        for index in range(len(lines)):
            
            # chi^2
            if 'chi^2:' in lines[index]:
                chi = lines[index].split()
                chi_square = float(chi[2])
            
            # Cosmology 
            if 'omega' in lines[index].split():
                cosmol = lines[index].split()
                h = float(cosmol[12])
                print(cosmol, h)
            
            # Lens
            if 'lens' in lines[index].split():
                lens= lines[index].split()
                vel_disp = float(lens[2])
                xl1 = float(lens[3])
                xl2 = float(lens[4])
                ql = float(lens[5])
                pa = float(lens[6])
                r_core = float(lens[7])


resdir.append({
    "losID" : system_id,
    "vel_disp" : mass_scale,
    "xl2" : xl1,
    "xl2" : xl2,
    "pa" : pa,
    "ql" : ql,
    "ys1" : ys1, 
    "ys2" : ys2,
    "chi_total" : chi_total,
    "h" : h,
    })


with open("./quasars_%s_nimgs_%s.csv" % (args["los"], args["nimgs"]), 'w') as fout:
    json.dump(resdir, fout)
