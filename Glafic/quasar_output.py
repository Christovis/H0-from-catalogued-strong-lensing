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
args["opt_explor"] = int(sys.argv[5])
outdir = os.fsencode(args["outdirstr"])

resdir = []  # initialize output dictionary

# Sort system_id
if args["opt_explor"] > 0:
    files_explore = []
    file_names = args["outdirstr"]+"fitH0_*_explore.dat"
    for file in glob.glob(file_names):
        files_explore.append(file)
    files_explore.sort(key=lambda ff: int(ff.split("_")[-2]))

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

files_ein = []
file_names = args["outdirstr"]+"fitH0_*_ein.dat"
for file in glob.glob(file_names):
    files_ein.append(file)
files_ein.sort(key=lambda ff: int(ff.split("_")[-2]))

#files_h = []
#file_names = args["outdirstr"]+"fitH0_*_vary.dat"
#for file in glob.glob(file_names):
#    files_h.append(file)
#files_h.sort(key=lambda ff: int(ff.split("_")[-2]))

# Run through files
for ff in range(num_of_files):
    system_id = int(files_point[ff].split('_')[-2])
    try:
        if args["opt_explor"] > 0:
            exp_indx = [ii for ii,fil in enumerate(files_explore) if int(fil.split('_')[-2]) == system_id][0]
        opt_indx = [ii for ii,fil in enumerate(files_optresult) if int(fil.split('_')[-2]) == system_id][0]
        ein_indx = [ii for ii,fil in enumerate(files_ein) if int(fil.split('_')[-2]) == system_id][0]
        #h_indx = [ii for ii,fil in enumerate(files_h) if int(fil.split('_')[-2]) == system_id][0]
    except:
        print('Failed to match', system_id, ff, files_explore[system_id,].split('_')[-2])
        continue
    
    # Point
    with open(files_point[ff], "r") as f:
        lines = f.readlines()

        if (int(args["nimgs"]) < len(lines)-2 or len(lines)-1 < int(args["nimgs"])):
            print("wrong nr. of images: losID %d; %d; %s" % (system_id, len(lines)-1, args["nimgs"]))
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

    # chi^2
    if args["opt_explor"] > 0:
        with open(files_explore[exp_indx], "r") as f:
            lines = f.readlines()
            try:
                indx = [ii for ii, ll in enumerate(lines) if "best-fit" in ll][0]
            except:
                print("Model parameter exceeded range")
                continue
            chi_square = (lines[indx].split()[-1])
            if float(chi_square) > 1e3:
                failed = True
            else:
                failed = False
        if failed is True:
            print("modelling failed for system %d" % system_id)
            continue
    
    # Optimization
    with open(files_optresult[opt_indx], "r") as f:
        lines = f.readlines()
        if args["opt_explor"] > 0:
            indx = [ii for ii, ll in enumerate(lines) if ('chi^2' in ll)][0]
        elif args["opt_explor"] == 0:
            chi_squares = np.asarray([[float(ll.split()[2]), ii] for ii, ll in enumerate(lines) if "chi^2" in ll])
            chi_squares = chi_squares.T
            indx = int(chi_squares[1, np.argmin(chi_squares[0,:])])
            chi_square = np.argmin(chi_squares[0, :])
        
        # Cosmology 
        if 'hubble' in lines[indx+6].split():
            cosmol = lines[indx+6].split()
            hubble = float(cosmol[-1])
        
        # Lens
        if 'lens   sie' in lines[indx+8]:
            lens= lines[indx+8].split()
            vel_disp = float(lens[2])
            xl1 = float(lens[3])
            xl2 = float(lens[4])
            ql = float(lens[5])
            pa = float(lens[6])
            r_core = float(lens[7])
            
        if args["los"] == "with_los":
            # external perturbation
            if 'lens   pert' in lines[indx+9]:
                pert = lines[indx+9].split()
                shear = float(pert[5])
                shear_angle = float(pert[6])
        
            # Source
            if 'point  ' in lines[indx+10]:
                lens= lines[indx+10].split()
                ys1 = float(lens[2])
                ys2 = float(lens[3])
        elif args["los"] == "no_los":
            # Source
            if 'point  ' in lines[indx+9]:
                lens= lines[indx+9].split()
                ys1 = float(lens[2])
                ys2 = float(lens[3])
    
    # Einstein-Ring
    with open(files_ein[ein_indx], "r") as f:
        lines = f.readlines()
        if not lines: 
            print("no Einstein radius for system %d" % system_id)
            continue
        Re = float(lines[0].split()[2])

    ## Hubble constant
    #with open(files_h[h_indx], "r") as f:
    #    lines = f.readlines()
    #    chi_squares = np.asarray([[float(ll.split()[1]), ii] for ii, ll in enumerate(lines)])
    #    chi_squares = chi_squares.T
    #    indx = int(chi_squares[1, np.argmin(chi_squares[0,:])])
    #    chi_square = np.argmin(chi_squares[0, :])
    #    hubble = float(lines[indx].split()[0])
    
    if args["los"] == "no_los":
        if args["nimgs"] == "2":
            resdir.append({
                "losID" : system_id,
                "vel_disp" : vel_disp,
                "xl2" : xl1,
                "xl2" : xl2,
                "pa" : pa,
                "ql" : ql,
                "r_core" : r_core,
                "ys1" : ys1, 
                "ys2" : ys2,
                "chi_total" : float(chi_square),
                "Re" : Re,
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
        elif args["nimgs"] == "4":
            resdir.append({
                "losID" : system_id,
                "vel_disp" : vel_disp,
                "xl2" : xl1,
                "xl2" : xl2,
                "pa" : pa,
                "ql" : ql,
                "r_core" : r_core,
                "ys1" : ys1, 
                "ys2" : ys2,
                "chi_total" : float(chi_square),
                "Re" : Re,
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
    elif args["los"] == "with_los":
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
                "chi_total" : float(chi_square),
                "Re" : Re,
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
        elif args["nimgs"] == "4":
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
                "chi_total" : float(chi_square),
                "Re" : Re,
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


with open("./quasars_%s_nimgs_%s_%s.json" % (args["los"], args["nimgs"], args["version"]), 'w') as fout:
    json.dump(resdir, fout)

