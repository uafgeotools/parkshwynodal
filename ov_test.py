import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

jet = ['B737', 'B738', 'B739']
Turboprop = ['DH8A','B190','BE20','C208','PC12','DH3T']
piston = ['C185','C182','C206','DHC2','GA8','PA31']
Heli = ['R44']

# Move to suplumentary material
piston = ['CH7B', 'PA30', 'PA32', 'C172','C180']
Turboprop = ['B18T','C441','AT73','SW4']
jet = ['B733', 'B763', 'B772', 'B77W', 'B788', 'B789']

eq = 'R44'

# Define the directory where your files are located
file = 'output/inv_results_old/' + eq + '_full_inv_results.csv' 

with open(file, 'r') as f:
    # Read the data from the file and append it to the list
    data = []
    fdiff = []
    for line in f.readlines():
        lines = line.split(',')

        if lines[-2] == "Forward Model" or lines[-5] == "00":
            continue

        peaks = np.array(lines[9])
        peaks = str(peaks) 
        peaks = peaks.replace('[', '').replace(']', '')
        peaks = np.array(peaks.split(' '))
        peak_old = 0
        for peak in peaks:
            if peak == '':
                continue
            peak = float(peak)
            if len(peaks) == 0 or peak == peaks[0]:
                peak_old = float(peak)
                continue

            diff = float(peak) - float(peak_old)

            fdiff.append(diff)
            peak_old = float(peak)
            data.append(peak)

med = np.median(fdiff)
if eq == 'C182':
    med = 18
elif eq == 'C185':
    med = 19.5
elif eq == 'DH8A':
    med = 15
elif eq == 'BE20':
    med = 27.5
    #med = 25.5
elif eq == 'C208':
    med = 29
elif eq == 'CH7B':
    med = 15
elif eq == 'B190':
    med = 24.5
elif eq == 'PC12':
    med = 28
elif eq == 'DH3T':
    med = 13
    med = 26
elif eq == 'C206':
    med = 19
elif eq == 'DHC2':
    med = 17.5
elif eq == 'GA8':
    med = 20
elif eq == 'PA31':
    med = 18.5
elif eq == 'R44':
    med = 13.35
elif eq == 'B739':
    med = 35.5
elif eq == 'B738':
    med = 73
elif eq == 'B737':
    med = 68
plt.figure(figsize=(10, 5))
bins = np.arange(min(data), max(data) + 3, 3)
plt.hist(data, color='k', bins=bins, alpha=0.5, edgecolor='black')
plt.text(0.99, 0.95, eq, transform=plt.gca().transAxes, fontsize=10, va='top', ha='right')
plt.xlim(0,275)
for g in range(0,50):
    plt.axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)

plt.show()
