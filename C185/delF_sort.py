import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = open('C185data_atmosphere.txt', 'r')
file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

plt.figure()

# Create a dictionary to store the color for each tail number
color_dict = {}
y_pos_dict = {}

# Create a dictionary to store the peaks for each tail number
med_dict = {}
date_dict = {}
count = 0
peak_old = None 
med_store = []
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]
    count += 1
    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    f1 = []
    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        if peak_old is None:
            peak_old = peak  # Set peak_old for the first iteration
            continue
        diff = float(peak) - float(peak_old)
        #if diff > 21 or diff < 18:
        #    continue
        f1.append(diff)
        peak_old = peak

    if not np.isnan(np.nanmedian(f1)):
        all_med = np.nanmedian(f1)
        med_store.append(all_med)
    print(flight_num)
    print(all_med)
    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            if str(tail_num) != '10512184':
                continue
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)
                y_pos_dict[tail_num] = []
                med_dict[tail_num] = []
                date_dict[tail_num] = []
            y_pos_dict[tail_num].extend([count])
            med_dict[tail_num].extend([all_med])
            date_dict[tail_num].extend([lines[3]])
for tail_num, med in med_dict.items():

    color = color_dict[tail_num]
    dates = date_dict[tail_num]
    y =  y_pos_dict[tail_num]   
    plt.scatter(med, y, c=color,label=tail_num)
plt.legend()
plt.xlim(0,80)
plt.show()

plt.figure()
for tail_num, med in med_dict.items():
    color = color_dict[tail_num]
    dates = date_dict[tail_num]
    y =  y_pos_dict[tail_num] 
    if tail_num == 10572742:  
        plt.hist(med, bins=1000, color=color, label=tail_num)
    elif tail_num == 10512184:
        plt.hist(med, bins=800, color=color, label=tail_num)
    else:
        plt.hist(med, bins=3, color=color, label=tail_num)
plt.legend()
plt.xlim(0,80)
plt.show()

plt.figure()
plt.hist(med_store, bins=2000)
plt.xlim(0,80)
plt.show()
