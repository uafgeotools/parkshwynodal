import numpy as np
import pandas as pd
from datetime import datetime, timezone
from prelude import *

def add_arrival_time(file_path, arrival_time):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            line = line.strip()
            if line:
                file.write(line + str(arrival_time) + ',\n')
        file.close()

equip = 'C185'
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

sta_f = open('input/all_station_crossing_db_' + equip + '.txt','r')
second_column = []
for line in sta_f.readlines():
    val = line.split(',')
    if len(val) >= 2:
        second_column.append(val[1])
sta_f.close()
second_column_array = np.array(second_column)

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')
for li in file_in.readlines():
    text = li.split(',')
    flight_num = text[1]
    if flight_num not in second_column_array:
        continue
    date = text[0]
    try:
        sta = int(text[9])
    except:
        continue

    c = 343

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
 
    timestamps = flight_data['snapshot_id']
    speed = flight_data['speed']
    altitude = flight_data['altitude']

    closest_x, closest_y, dist_km, closest_time, tarrive, alt, sp, elevation, speed_mps, height_m, dist_m, tmid = closest_approach_UTM(seismo_latitudes, seismo_longitudes, flight_latitudes, flight_longitudes, timestamps, altitude, speed, stations, elevations, c, sta)

    ta_old = calc_time(tmid,dist_m,height_m,343)
    
    ht = datetime.fromtimestamp(ta_old, tz=timezone.utc)
    month = ht.month
    day = ht.day		

    if len(str(day)) == 1:
        day = '0'+str(day)

    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
    if Path(file_name).exists():
        add_arrival_time(file_name, ta_old)
    else:
        print(file_name)
        continue

    output2 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/overtonepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
    if Path(output2).exists():
       add_arrival_time(output2, ta_old)
    else:
        print(file_name)
        continue

    output3 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/timepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
    if Path(output3).exists():
        add_arrival_time(output3, ta_old)

    else:
        print(file_name)
        continue