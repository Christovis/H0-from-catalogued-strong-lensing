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

    if filename.endswith(".best"):
        print(filename)
        system_id = filename.split('_')[-1].split('.')[0]
        with open(outdirstr+filename, "r") as f:
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
                    ys1 = float(source[1])
                    ys2 = float(source[2])
                
                if 'CHISQ:' in lines[index].split():
                    chi_square = lines[index].split()
                    chi_total = float(chi_square[1])

                elif 'h:' in lines[index].split():
                    h = float(lines[index].split()[1])

                elif 'images:' in lines[index].split():
                    
                    if args["nimgs"] == "2":
                        image_1 = lines[index+1].split()
                        x1_1 = float(image_1[9]) 
                        x2_1 = float(image_1[10]) 
                        mu_1 = float(image_1[11]) 
                        image_2 = lines[index+2].split()
                        x1_2 = float(image_2[9]) 
                        x2_2 = float(image_2[10]) 
                        mu_2 = float(image_2[11]) 
                        print('image 1 pos:', x1_1, x2_1)
                        print('image 2 pos:', x1_2, x2_2)
                    
                    elif args["nimgs"] == "4":
                        image_1 = lines[index+1].split()
                        x1_1 = float(image_1[9]) 
                        x2_1 = float(image_1[10]) 
                        mu_1 = float(image_1[11]) 
                        image_2 = lines[index+2].split()
                        x1_2 = float(image_2[9]) 
                        x2_2 = float(image_2[10]) 
                        mu_2 = float(image_2[11]) 
                        image_3 = lines[index+3].split()
                        x1_3 = float(image_3[9]) 
                        x2_3 = float(image_3[10]) 
                        mu_3 = float(image_3[11]) 
                        image_4 = lines[index+4].split()
                        x1_4 = float(image_4[9]) 
                        x2_4 = float(image_4[10]) 
                        mu_4 = float(image_4[11]) 
        
        filename = "Rein_" + system_id
        with open(outdirstr+filename, "r") as f:
            lines = f.readlines()
            Re = lines[0].split()[0]

    else:
        continue
    
    if args["nimgs"] == "2":
        resdir.append({
            "mass_scale" : mass_scale,
            "xl2" : xl1,
            "xl2" : xl2,
            "pa" : pa,
            "ql" : ql,
            "shear" : shear,
            "shear_angle" : shear_angle,
            "ys1" : ys1, 
            "ys2" : ys2,
            "chi_total" : chi_total,
            "Re" : Re,
            "h" : h,
            "x1_1" : x1_1, 
            "x2_1" : x2_1, 
            "x1_2" : x1_2,
            "x2_2" : x2_2,
            "mu_1" : mu_1, 
            "mu_2" : mu_2,
            })
    elif args["nimgs"] == "4":
        resdir.append({
            "mass_scale" : mass_scale,
            "xl2" : xl1,
            "xl2" : xl2,
            "pa" : pa,
            "ql" : ql,
            "shear" : shear,
            "shear_angle" : shear_angle,
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
            "mu_1" : mu_1, 
            "mu_2" : mu_2,
            "mu_3" : mu_3, 
            "mu_4" : mu_4,
            })


with open("./quasars_%s_%s.csv" % (args["los"], args["nimgs"]), 'w') as fout:
    json.dump(resdir, fout)

#lens_pos = np.asarray(lens_pos).T
#print(np.median(hubble_const))
#print(len(mass_scale), len(lens_pos[0, :]), len(hubble_const),
#      np.shape(ellipticity), np.shape(position_angle), np.shape(shear_angle))
#df = pd.DataFrame({
#    'MS' : mass_scale,
#    'LPX' : lens_pos[0, :],
#    'LPY' : lens_pos[1, :],
#    'ELL' : ellipticity,
#    'PA' : position_angle,
#    'SHE' : shear,
#    'SA' : shear_angle,
#    'H0' : hubble_const,
#})
#df.to_csv("./quasars_with_los.csv", index=True)

