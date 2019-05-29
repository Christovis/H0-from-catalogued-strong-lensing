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
args["version"] = sys.argv[4]
outdir = os.fsencode(args["outdirstr"])

resdir = []  # initialize output dictionary

# Run through files
for file in os.listdir(outdir):
    filename = os.fsdecode(file)

    if filename.endswith(".best"):
        system_id = filename.split('_')[-1].split('.')[0]
        with open(args["outdirstr"]+filename, "r") as f:
            lines = f.readlines()
    
            # Run through lines in file
            for index in range(len(lines)):
                
                if 'LENS' in lines[index].split():
                    lens = lines[index+1].split()
                    mass_scale = float(lens[1])
                    xl1 = float(lens[2])
                    xl2 = float(lens[3])
                    ql = float(lens[4])
                    pa = float(lens[5])
                    shear = float(lens[6])
                    shear_angle = float(lens[7])
                
                if 'SOURCE' in lines[index].split():
                    source = lines[index+1].split()
                    fluxs = float(source[1])
                    ys1 = float(source[2])
                    ys2 = float(source[3])
                
                if 'CHISQ:' in lines[index].split():
                    chi_square = lines[index].split()
                    chi_total = float(chi_square[1])

                elif 'h:' in lines[index].split():
                    h = float(lines[index].split()[1])

                
                elif 'images:' in lines[index].split():
                    if args["nimgs"] == "2":
                        # collect image 
                        image_a = lines[index+1].split()
                        image_b = lines[index+2].split()
                        x1 = [float(image_a[9]), float(image_b[9])]
                        x2 = [float(image_a[10]), float(image_b[10])]
                        mu = [float(image_a[11]), float(image_b[11])]
                        dt = [float(image_a[12]), float(image_b[12])]
                        
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
                        # collect image 
                        image_a = lines[index+1].split()
                        image_b = lines[index+2].split()
                        image_c = lines[index+3].split()
                        image_d = lines[index+4].split()
                        x1 = [float(image_a[9]), float(image_b[9]),
                              float(image_c[9]), float(image_d[9])]
                        x2 = [float(image_a[10]), float(image_b[10]),
                              float(image_c[10]), float(image_d[10])]
                        mu = [float(image_a[11]), float(image_b[11]),
                              float(image_c[11]), float(image_d[11])]
                        dt = [float(image_a[12]), float(image_b[12]),
                              float(image_c[12]), float(image_d[12])]
                        
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
                        x2_4 = x1[index[3]]
                        mu_4 = mu[index[3]]
                        dt_4 = dt[index[3]]

        filename = "Rein_" + system_id
        with open(args["outdirstr"]+filename, "r") as f:
            lines = f.readlines()
            Re = lines[0].split()[0]
        
        """
        filename = "magnification_" + system_id + ".dat"
        with open(args["outdirstr"]+filename, "r") as f:
            lines = f.readlines()
            
            # If system not correct number of multiple images
            if (len(lines)-2) < int(args["nimgs"]):
                continue

            if args["nimgs"] == "2":
                x1_1 = lines[2].split()[0]
                x2_1 = lines[2].split()[1]
                mu_1 = lines[2].split()[2]
                
                x1_2 = lines[3].split()[0]
                x2_2 = lines[3].split()[1]
                mu_2 = lines[3].split()[2]

            elif args["nimgs"] == "4":
                x1_1 = lines[2].split()[0]
                x2_1 = lines[2].split()[1]
                mu_1 = lines[2].split()[2]
                
                x1_2 = lines[3].split()[0]
                x2_2 = lines[3].split()[1]
                mu_2 = lines[3].split()[2]
                
                x1_3 = lines[4].split()[0]
                x2_3 = lines[4].split()[1]
                mu_3 = lines[4].split()[2]
                
                x1_4 = lines[5].split()[0]
                x2_4 = lines[5].split()[1]
                mu_4 = lines[5].split()[2]
        """
    else:
        continue
    
    if args["nimgs"] == "2":
        resdir.append({
            "losID" : system_id,
            "mass_scale" : mass_scale,
            "xl2" : xl1,
            "xl2" : xl2,
            "pa" : pa,
            "ql" : ql,
            "shear" : shear,
            "shear_angle" : shear_angle,
            "flux_src" : fluxs, 
            "ys1" : ys1, 
            "ys2" : ys2,
            "chi_total" : chi_total,
            "Re" : Re,
            "h" : h,
            "x1_1" : x1_1, 
            "x2_1" : x2_1, 
            "x1_2" : x1_2,
            "x2_2" : x2_2,
            #"mags_1" : mags_1, 
            #"mags_2" : mags_2,
            "mu_1" : mu_1, 
            "mu_2" : mu_2,
            "dt_1" : dt_1, 
            "dt_2" : dt_2,
            })
    elif args["nimgs"] == "4":
        resdir.append({
            "losID" : system_id,
            "mass_scale" : mass_scale,
            "xl2" : xl1,
            "xl2" : xl2,
            "pa" : pa,
            "ql" : ql,
            "shear" : shear,
            "shear_angle" : shear_angle,
            "flux_src" : fluxs, 
            "ys1" : ys1, 
            "ys2" : ys2,
            "chi_total" : chi_total,
            "Re" : Re,
            "h" : h,
            "x1_1" : x1_1, 
            "x2_1" : x2_1, 
            "x1_2" : x1_2,
            "x2_2" : x2_2,
            "x1_3" : x1_3, 
            "x2_3" : x2_3, 
            "x1_4" : x1_4,
            "x2_4" : x2_4,
            #"mags_1" : mags_1, 
            #"mags_2" : mags_2,
            #"mags_3" : mags_3, 
            #"mags_4" : mags_4,
            "mu_1" : mu_1, 
            "mu_2" : mu_2,
            "mu_3" : mu_3, 
            "mu_4" : mu_4,
            "dt_1" : dt_1, 
            "dt_2" : dt_2,
            "dt_3" : dt_3, 
            "dt_4" : dt_4,
            })


with open("./quasars_%s_nimgs_%s_%s.json" % (args["los"], args["nimgs"], args["version"]), 'w') as fout:
    json.dump(resdir, fout)

