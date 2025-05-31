import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = open('C185data_atm_full.txt', 'r')
file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']
alts = file2['ELEVATION']
time = file2['TIME']

# Create a dictionary to store the color for each tail number
color_dict = {}
y_pos_dict = {}

# Create a dictionary to store the peaks for each tail number
peaks_dict = {}
date_dict = {}

count = 0
xx = 'yes'
all_med = {}
y_med = {}
mad = {}
date_med = {}
flig = {}
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]
    nodes = int(lines[2])
    distance = float(lines[6])
    vo = float(lines[4]) #try with vo and also with diffrence between heading direction and travel direction
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/'+str(lines[0]) + '_positions/' + str(lines[0]) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    t = flight_data['snapshot_id']
    speed = flight_data['speed']
    #alt = flight_data['altitude']

    #time_diff = t - float(lines[3])  # Subtract lines[3] from each value in the t array

    # Subtract the value of lines[3] from each value in the time_diff array
    time_diff = np.subtract(t, float(lines[3]))

    closest_time_index = np.argmin(np.abs(time_diff))
    closest_time = t[closest_time_index]
    #closest_alt = alt[closest_time_index]
    closest_speed = speed[closest_time_index]
    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    count += 1

    ppp = []
    date = []
    y=[]
    f1 = []
    peak_old = 0
    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            if np.abs(float(peak) - float(peak_old))< 10:
                continue
            ppp.append(float(peak))
            #date.append(float(lines[3]))
            date.append(lines[3])
            y.append(count)
            if len(peaks) == 0 or peak == peaks[0]:
                peak_old = peak
                continue
            diff = float(peak) - float(peak_old)
            if diff > 22 or diff < 18:
                continue
            f1.append(diff)
        peak_old = float(peak)

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            #tail_num = alts[lp] * 0.3048  # Convert altitude from feet to meters
            tail_num = flight_num #tail_nums[lp]
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = [] #alts[lp] * 0.3048 # np.random.rand(3,)

                peaks_dict[tail_num] = []
                y_pos_dict[tail_num] = []
                date_dict[tail_num] = []

                all_med[tail_num] = []
                y_med[tail_num] = []
                mad[tail_num] = []
                date_med[tail_num] = []
                flig[tail_num] = []
            for i in range(len(ppp)):
                color_dict[tail_num].extend([closest_speed]) #[alts[lp] * 0.3048])
            peaks_dict[tail_num].extend(ppp)
            date_dict[tail_num].extend(date)
            y_pos_dict[tail_num].extend(y)

            all_med[tail_num].extend([np.nanmedian(f1)])
            #date_med[tail_num].extend([float(lines[3])])
            date_med[tail_num].extend([lines[3]])
            y_med[tail_num].extend([count])
            mad[tail_num].extend([np.median(np.absolute(f1 - np.median(f1)))])
    flig[tail_num].extend([flight_num])
fig,ax1 = plt.subplots(1, 1, sharex=False, figsize=(8,6))     

ax1.margins(x=0)
#ax1.grid(axis='both') 
ax2 = fig.add_axes([0.83, 0.11, 0.07, 0.77], sharey=ax1)
ax3 = fig.add_axes([0.90, 0.11, 0.07, 0.77], sharey=ax1) 
ax1.set_title('Frequency Peaks')
# Sort tail numbers by their first color value in ascending order
sorted_tail_nums = sorted(color_dict.keys(), key=lambda tn: color_dict[tn][0] if color_dict[tn] else float('inf'))

for tail_num, peaks in peaks_dict.items():
    print(str(flig[tail_num]))
    #if str(flight_num) == '529754214' or str(tail_num) == '528698927' or str(tail_num) == '529409728' or str(tail_num) == '529416700':
    #    go ='y'
 
    #if '528698927' in flig[tail_num]:
    #    #go ='y'
    #else:
    #    continue
    print(peaks_dict[tail_num])
    color = color_dict[tail_num]
    dates = date_dict[tail_num]
    y =  y_pos_dict[tail_num]
    ax1.scatter(peaks, y, c=color) #, cmap='pink', edgecolors='black') #label=tail_num)
for tail_num, med in all_med.items():
    #if str(tail_num) == '529754214' or str(tail_num) == '528698927' or str(tail_num) == '529409728' or str(tail_num) == '529416700':
    #    go ='y'
 
    color = color_dict[tail_num][0:len(med)]  # Use the first color for the median points
    y= y_med[tail_num]
    dates = date_med[tail_num]
    ax2.scatter(med, y, c=color) #, cmap='pink', edgecolors='black') 
    m = mad[tail_num]
    ax3.scatter(m, y, c=color) #, cmap='pink', edgecolors='black')  
ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)

ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax2.set_title('Median')
ax3.set_title('MAD')
ax1.grid(axis='both') 
ax2.grid(axis='both') 
ax3.grid(axis='both') 
ax1.set_xlabel('Frequency')
ax2.set_xlabel('\u0394'+'F')
ax1.legend(loc='upper left',fontsize = 'x-small')
ax1.set_xlim(0, 300)
ax1.set_xticks(range(0, 251, 10)) 

#ax1.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
plt.show()
