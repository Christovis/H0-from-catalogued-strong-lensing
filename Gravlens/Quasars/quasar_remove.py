from __future__ import division
import os
import numpy as np
import pandas as pd

los="no_los"
directoryname = "/cosma5/data/dp004/dc-beck3/StrongLensing/Gravlens/" + \
                "Quasars/results/"+los+"/" 
directory = os.fsencode(directoryname)  # only python 3

# Run through files
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".best"):
        text_file = open(directoryname+filename, "r")
        
        # Run through lines in file
        for line in text_file:
            words = line.split(' ')
            if 'h:' in words:
                if np.abs(float(words[-1]) - 0.7) > 0.05:
                    print(filename)
                    print('H0', float(words[-1]), filename[:-4])
                    #os.remove(directoryname+filename)
                    #os.remove(directoryname+filename[:-4]+"dat")
                    #os.remove(directoryname+filename[:-4]+"startN")
                    #os.remove(directoryname+filename[:-4]+"chi")
                    #os.remove(directoryname+filename[:-4]+"start")

    else:
        continue

