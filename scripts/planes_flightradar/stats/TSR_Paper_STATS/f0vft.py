import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj
from matplotlib.gridspec import GridSpec
from matplotlib import colors as mcolors
import os
from doppler_funcs import calc_ft
file_list =[]
# loop through the directories in the directory
for file in os.listdir('/home/irseppi/REPOSITORIES/parkshwynodal/output/Inversion_Results/'):
    dir_path = os.path.join('/home/irseppi/REPOSITORIES/parkshwynodal/output/Inversion_Results/', file)
    if os.path.isfile(dir_path):
        file_list.append(dir_path)
diff = []
# Process each file in the list
for gg, file in enumerate(file_list):
    if os.path.getsize(file) == 0:
        continue
    file = open(file, 'r')
    # Parse the file line by line
    for line in file.readlines():
        lines = line.split(',')
        tprime0 = float(lines[3])
        v0 = float(lines[5])
        l = float(lines[6])
        c = float(lines[10])
        #if c >= v0:
        #    continue
        peaks_new = np.array(lines[7])
        peaks_new = np.array(lines[7])
        peaks_new = str(peaks_new)
        peaks_new = np.char.replace(peaks_new, '[', '')
        peaks_new = np.char.replace(peaks_new, ']', '')
        peaks_new = str(peaks_new)
        peaks_new = np.array(peaks_new.split(' '))
        for i in range(len(peaks_new)):
            try:
                f0 = float(peaks_new[i])
                tprime = tprime0
                t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                ft = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
                diff.append(ft - f0)
            except ValueError:
                continue

plt.figure(figsize=(10, 6))
plt.hist(diff, bins=100, color='blue', range=(0, 100))
median_val = np.nanmedian(diff)
plt.axvline(median_val, color='red', linestyle='dotted', linewidth=2, label=f'Median: {median_val:.2f}')
print(f"Median value: {median_val}")

mean_val = np.nanmean(diff)
plt.axvline(mean_val, color='green', linestyle='dashed', linewidth=2, label=f'Mean: {mean_val:.2f}')
print(f"Mean value: {mean_val}")

plt.legend()
plt.show()
