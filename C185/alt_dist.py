import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj

file = open('output3.txt', 'r')
file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.txt', sep=",")
tail_nums = file2['10512184']
flights = file2['527958214']
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

list_f = []
color_dict = {}
label_dict = {}

MODE = 1
mode1_flight = []
mode2_flight = []

plt.figure()
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = int(lines[1])
    if flight_num in list_f:
        continue
    else:
        list_f.append(flight_num)

    peaks = np.array(lines[7])
    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    for peak in peaks:
        try:
            peak = float(peak[0:-1])
        except:
            continue
        if peak <111 or peak > 126:
            continue
        else:
            if abs(float(peak) - 124.5) <= abs(float(peak) - 114):
                mode1_flight.append(flight_num)
                break
            elif abs(float(peak) - 114) <= abs(float(peak) - 124.5):
                mode2_flight.append(flight_num)
                break

    for lp in range(len(flights)):
        if flight_num == int(flights[lp]):
            tail_num = tail_nums[lp]
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)

    date = lines[0]
 
    month = date[4:6]
    day = date[6:8]

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+ '_positions/2019'+month+day+ '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")

    altitude = flight_data['altitude']
    lat = flight_data['latitude']
    lon = flight_data['longitude']

    alts = []
    dists = []

    for a in range(1,len(altitude)-1):
        if a == 1:
            dist  = 0
            x1,y1 = utm_proj(lon[a], lat[a])
        else:
            x,y = utm_proj(lon[a], lat[a])
            dist = np.sqrt((x-x1)**2 + (y-y1)**2) + dist
            x1 = x
            y1 = y

        alts.append(altitude[a])
        dists.append(dist/1000)

    if MODE == 1:
        print(mode1_flight)
        title = 'Mode 1'
        if flight_num in mode1_flight:
            if tail_num not in label_dict:
                plt.plot(dists, alts, color=color_dict[tail_num], label=tail_num)
                label_dict[tail_num] = tail_num
            else:
                plt.plot(dists, alts, color=color_dict[tail_num])
        else:
            continue
    elif MODE == 2:
        print(mode2_flight)
        title = 'Mode 2'
        if flight_num in mode2_flight:
            if tail_num not in label_dict:
                plt.plot(dists, alts, color=color_dict[tail_num], label=tail_num)
                label_dict[tail_num] = tail_num
            else:
                plt.plot(dists, alts, color=color_dict[tail_num])
        else:
            continue
plt.legend() 
plt.title(title)
plt.show()
