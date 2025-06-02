import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = open('C185data_atm_full.txt', 'r')

file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

# Create a dictionary to store the color for each tail number
color_dict = {}
y_pos_dict = {}
peaks_dict = {}
all_med = {}
mad = {}

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

    count += 1

    ppp = []
    y=[]
    f1 = []
    peak_old = 0
    for peak in peaks:
        if np.abs(float(peak) - float(peak_old))< 10:
            continue

        ppp.append(float(peak))
        y.append(count)

        if len(peaks) == 0 or peak == peaks[0]:
            peak_old = float(peak)
            continue

        diff = float(peak) - float(peak_old)
        if diff > 21 or diff < 18:
            continue
        f1.append(diff)
        peak_old = float(peak)

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,) 
                peaks_dict[tail_num] = []
                y_pos_dict[tail_num] = []
                all_med[tail_num] = []
                mad[tail_num] = []
                break
        else:
            continue
    peaks_dict[tail_num].extend(ppp)
    y_pos_dict[tail_num].extend(y)
    all_med[tail_num].extend([np.nanmedian(f1)])
    mad[tail_num].extend([np.median(np.absolute(f1 - np.median(f1)))])

fig,ax1 = plt.subplots(1, 1, sharex=False, figsize=(8,6))     

ax1.margins(x=0)
ax2 = fig.add_axes([0.83, 0.11, 0.07, 0.77], sharey=ax1)

pos = 1
tail_num_hold = 0
color_dict[10572742] = [1.0, 0.5, 0.0]  # Orange color in RGB
color_dict[10512184] = [0.0, 0.5, 1.0]  # Blue color in RGB
for tail_num, peaks in peaks_dict.items():
    color = color_dict[tail_num]
    if str(tail_num) != '10572742' and str(tail_num) != '10512184':
        continue

    new_y = []
    y_m = [pos]
    for i,index in enumerate(y_pos_dict[tail_num]):
        if i == 0:
            index_old = index
        elif index == index_old:
            index_old = index
        else:
            index_old = index
            pos += 1
            y_m.append(pos)
        new_y.append(pos)
        if i == len(y_pos_dict[tail_num]) - 1:
            pos += 1
    med = all_med[tail_num]
    md = mad[tail_num]
    #ax1.scatter(peaks, new_y, c=color, label=tail_num)  
    ax1.hist(peaks, bins=200, color=color, alpha=0.8, label=tail_num)  # Use the first color for the histogram

    #ax2.scatter(med, y_m, c=color)
    ax2.hist(med, bins=10, color=color, alpha=0.8)  # Use the first color for the histogram


ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax1.grid(which='both', axis='both', ls='--') 
ax2.grid(axis='both', linestyle='--') 
ax1.set_xlabel('Frequency')
ax2.set_xlabel('Median '+'\u0394'+'F')
ax1.legend(loc='upper left',fontsize = 'large')
ax1.set_xlim(10, 300)
ax1.set_xticks(range(10, 270, 10)) 
ax2.set_yticks(range(0, 100, 10))

plt.show()
