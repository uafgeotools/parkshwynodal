import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = open('C185data_atm_full.txt', 'r')

file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

# Create a dictionary to store the color for each tail number
color_dict = {}
peaks_dict = {}
all_med = {}

count = 0
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]

    peaks = np.array(lines[7])
    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')
    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))


    ppp = []
    f1 = []
    peak_old = 0
    for peak in peaks:
        if np.abs(float(peak) - float(peak_old))< 10:
            continue
        ppp.append(float(peak))

        if len(peaks) == 0 or peak == peaks[0]:
            peak_old = float(peak)
            continue

        diff = float(peak) - float(peak_old)
        #if diff > 21 or diff < 18:
        #    continue
        f1.append(diff)
        peak_old = float(peak)

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = []
                peaks_dict[tail_num] = []
                all_med[tail_num] = []
                break
        else:
            continue
    peaks_dict[tail_num].extend(ppp)
    all_med[tail_num].extend([np.nanmedian(f1)])

fig,ax1 = plt.subplots(1, 1, sharex=False, figsize=(50,20))     

ax1.margins(x=0)
ax2 = fig.add_axes([0.83, 0.11, 0.1, 0.77], sharey=ax1)

pos = 1
tail_num_hold = 0
color_dict[10572742] = [1.0, 0.5, 0.0]  # Orange color in RGB
color_dict[10512184] = [0.0, 0.5, 1.0]  # Blue color in RGB

for tail_num, peaks in peaks_dict.items():
    color = color_dict[tail_num]
    med = all_med[tail_num]
    if str(tail_num) != '10572742' and str(tail_num) != '10512184':
        continue

    ax1.hist(peaks, bins=270, color=color, alpha=0.8, label=tail_num, zorder = 10)  
    ax2.hist(med, bins=270, color=color, alpha=0.8, zorder = 10)  
    ax1.hist(peaks, bins=270, color=color, histtype='step',zorder = 15)  
    ax2.hist(med, bins=270, color=color, histtype='step', zorder = 15)  

ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax1.set_xlabel('Frequency', fontsize=16)
ax2.set_xlabel('Median '+'\u0394'+'F', fontsize=16)
ax1.legend(loc='upper left',fontsize = 'x-large')
ax1.set_xlim(10, 298)
ax1.set_xticks(range(10, 280, 10)) 
ax1.set_yticks(range(0, 90, 10))
ax1.tick_params(axis='both', labelsize=14)  # Increase font size for tick labels
ax2.tick_params(axis='both', labelsize=14)  # Increase font size for tick labels
ax2.set_xticks(np.arange(18.5, 22, 1))
ax2.set_xlim(18, 22)
ax1.set_ylim(0, 80)
ax1.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
ax2.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
del_f_t1 = 19.62
del_f_t2_1 = 19.17
del_f_t2_2 = 20.56
for g in range(0,14):
    ax1.axvline(x= (1 + g) * del_f_t1, color = [1.0, 0.5, 0.0], ls = '--', zorder=0, linewidth=1)
    ax1.axvline(x= (1 + g) * del_f_t2_1, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
    ax1.axvline(x= (1 + g) * del_f_t2_2, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
ax2.axvline(x=del_f_t1, color = [1.0, 0.5, 0.0], ls = '--', zorder=0, linewidth=1)
ax2.axvline(x=del_f_t2_1, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
ax2.axvline(x=del_f_t2_2, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
plt.show()
