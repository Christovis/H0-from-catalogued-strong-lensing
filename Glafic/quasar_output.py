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
args["version"] = sys.argv[4]
outdir = os.fsencode(args["outdirstr"])

resdir = []  # initialize output dictionary

# Sort system_id
files_optresult = []
file_names = args["outdirstr"]+"fitH0_*_optresult.dat"
for file in glob.glob(file_names):
    files_optresult.append(file)
files_optresult.sort(key=lambda ff: int(ff.split("_")[-2]))

files_point = []
file_names = args["outdirstr"]+"fitH0_*_point.dat"
for file in glob.glob(file_names):
    files_point.append(file)
files_point.sort(key=lambda ff: int(ff.split("_")[-2]))

num_of_files = len(glob.glob(file_names))
# Run through files
for ff in range(num_of_files):
    system_id = ff

    with open(files_point[system_id], "r") as f:
        lines = f.readlines()
        print(files_point[system_id])
        print(lines)

        if len(lines) < int(args["nimgs"])+1:
            print("wrong nr. of images for system %d" % system_id)
            continue

        if args["nimgs"] == "2":
            image_a = lines[1].split()
            image_b = lines[2].split()
            x1 = [float(image_a[0]), float(image_b[0])]
            x2 = [float(image_a[1]), float(image_b[1])]
            mu = [float(image_a[2]), float(image_b[2])]
            dt = [float(image_a[3]), float(image_b[3])]

            # sort images by arrival time
            index = np.argsort(dt)

            # assign images
            x1_1 = x1[index[0]]
            x2_1 = x2[index[0]]
            mu_1 = mu[index[0]]
            dt_1 = dt[index[0]]
            x1_2 = x1[index[1]]
            x2_2 = x2[index[1]]
            mu_2 = mu[index[1]]
            dt_2 = dt[index[1]]
            
        elif args["nimgs"] == "4":
            image_a = lines[1].split()
            image_b = lines[2].split()
            image_c = lines[3].split()
            image_d = lines[4].split()
            print(image_d)
            x1 = [float(image_a[0]), float(image_b[0]),
                  float(image_c[0]), float(image_d[0])]
            x2 = [float(image_a[1]), float(image_b[1]),
                  float(image_c[1]), float(image_d[1])]
            mu = [float(image_a[2]), float(image_b[2]),
                  float(image_c[2]), float(image_d[2])]
            dt = [float(image_a[3]), float(image_b[3]),
                  float(image_c[3]), float(image_d[3])]

            # sort images by arrival time
            index = np.argsort(dt)

            # assign images
            x1_1 = x1[index[0]]
            x2_1 = x2[index[0]]
            mu_1 = mu[index[0]]
            dt_1 = dt[index[0]]
            x1_2 = x1[index[1]]
            x2_2 = x2[index[1]]
            mu_2 = mu[index[1]]
            dt_2 = dt[index[1]]
            x1_3 = x1[index[2]]
            x2_3 = x2[index[2]]
            mu_3 = mu[index[2]]
            dt_3 = dt[index[2]]
            x1_4 = x1[index[3]]
            x2_4 = x2[index[3]]
            mu_4 = mu[index[3]]
            dt_4 = dt[index[3]]

    with open(files_optresult[system_id], "r") as f:
        lines = f.readlines()
        # Run through lines in file bottom-to-top
        for index in range(len(lines)-1, 0, -1):
            # End of results
            if 'optimize' in lines[index]:
                break 

            # chi^2
            if 'chi^2' in lines[index]:
                chi = lines[index].split()
                chi_square = float(chi[2])
            
            # Cosmology 
            if 'omega' in lines[index].split():
                cosmol = lines[index].split()
                hubble = float(cosmol[-1])
            
            # Lens
            if 'lens   sie' in lines[index]:
                lens= lines[index].split()
                vel_disp = float(lens[2])
                xl1 = float(lens[3])
                xl2 = float(lens[4])
                ql = float(lens[5])
                pa = float(lens[6])
                r_core = float(lens[7])
            
            # external perturbation
            if 'lens   pert' in lines[index]:
                pert = lines[index].split()
                shear = float(pert[5])
                shear_angle = float(pert[6])
            
            # Source
            if 'point  ' in lines[index]:
                lens= lines[index].split()
                ys1 = float(lens[2])
                ys2 = float(lens[3])

    if args["nimgs"] == "2":
        resdir.append({
            "losID" : system_id,
            "vel_disp" : vel_disp,
            "xl2" : xl1,
            "xl2" : xl2,
            "pa" : pa,
            "ql" : ql,
            "shear" : shear,
            "shear_angle" : shear_angle,
            "r_core" : r_core,
            "ys1" : ys1, 
            "ys2" : ys2,
            "chi_total" : chi_square,
            "h" : hubble,
            "x1_1" : x1_1,
            "x2_1" : x2_1,
            "mu_1" : mu_1,
            "dt_1" : dt_1,
            "x1_2" : x1_2,
            "x2_2" : x2_2,
            "mu_2" : mu_2,
            "dt_2" : dt_2,
            })
    if args["nimgs"] == "4":
        resdir.append({
            "losID" : system_id,
            "vel_disp" : vel_disp,
            "xl2" : xl1,
            "xl2" : xl2,
            "pa" : pa,
            "ql" : ql,
            "shear" : shear,
            "shear_angle" : shear_angle,
            "r_core" : r_core,
            "ys1" : ys1, 
            "ys2" : ys2,
            "chi_total" : chi_square,
            "h" : hubble,
            "x1_1" : x1_1,
            "x2_1" : x2_1,
            "mu_1" : mu_1,
            "dt_1" : dt_1,
            "x1_2" : x1_2,
            "x2_2" : x2_2,
            "mu_2" : mu_2,
            "dt_2" : dt_2,
            "x1_3" : x1_3,
            "x2_3" : x2_3,
            "mu_3" : mu_3,
            "dt_3" : dt_3,
            "x1_4" : x1_4,
            "x2_4" : x2_4,
            "mu_4" : mu_4,
            "dt_4" : dt_4,
            })
    print(system_id)
    print(resdir[-1])
    break


with open("./quasars_%s_nimgs_%s_%s.json" % (args["los"], args["nimgs"], args["version"]), 'w') as fout:
    json.dump(resdir, fout)

