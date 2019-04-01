from __future__ import division
import os
import numpy as np
import pandas as pd

directoryname = "/cosma5/data/dp004/dc-beck3/StrongLensing/Gravlens/Quasars/results/with_los/" 
directory = os.fsencode(directoryname)

mass_scale = []
lens_pos = []
ellipticity = []
position_angle = []
shear = []
shear_angle = []
hubble_const = []

# Run through files
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".best"):
        text_file = open(directoryname+filename, "r")
        
        # Run through lines in file
        for line in text_file:
            words = line.split(' ')
            if 'alpha' in words:
                mass_scale.append(float(words[5]))
                lens_pos.append([float(words[6]), float(words[7])])
                ellipticity.append(float(words[8]))
                position_angle.append(float(words[9]))
                shear.append(float(words[10]))
                shear_angle.append(float(words[11]))
            elif 'h:' in words:
                hubble_const.append(float(words[-1]))
    else:
        continue

lens_pos = np.asarray(lens_pos).T
print(np.median(hubble_const))
#print(len(mass_scale), len(lens_pos[0, :]), len(hubble_const),
#      np.shape(ellipticity), np.shape(position_angle), np.shape(shear_angle))
df = pd.DataFrame({
    'MS' : mass_scale,
    'LPX' : lens_pos[0, :],
    'LPY' : lens_pos[1, :],
    'ELL' : ellipticity,
    'PA' : position_angle,
    'SHE' : shear,
    'SA' : shear_angle,
    'H0' : hubble_const,
    })
df.to_csv("./quasars_with_los.csv", index=True)
