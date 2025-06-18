import numpy as np
import pandas as pd
import json
import obspy
import datetime
from datetime import datetime, timezone
from pyproj import Proj
from prelude import *
from scipy.signal import find_peaks, spectrogram
from plot_func import *
from obspy.clients.nrl import NRL
import os
#add a way to get the correct time and save it with picks
#add a way to only pick 5 images for each aircraft type

nrl = NRL()

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

rerun_fig = False #Flag rerun the figures without saving the inversion results = True
equip_count_dict = {}
tailnumber_dict = {}
# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
for li in file_in.readlines():
    text = li.split(',')
    date = text[0]
    flight_num = text[1]
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate
    dist_m = float(text[4])   # Distance in meters
    closest_time = float(text[5])
    alt = float(text[6]) 
    speed_mps = float(text[7])  # Speed in meters per second
    sta = text[9]
    equip = text[10]

    if equip[0:3] == 'B73' or equip == 'nan':
        print(equip[0:3])
        continue
    #if equip[0:1] == 'B7' and dist_m < 1000:
    #    continue

    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    spec_dir = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
    #'/scratch/irseppi/nodal_data/plane_info/' + folder_spec +'/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'
    
    if os.path.exists(spec_dir) and rerun_fig == False:
        continue
    
    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + date + '_flights.csv', sep=",")
    flight = flight_data['flight_id']
    flight = flight.values.tolist()
    tailnumber = flight_data['aircraft_id']

    # get the index of the flight equivalent to the flight number
    index = flight.index(int(flight_num))
    #Fix this section to use files to count tailnumbers so you can get accurate counts
    if tailnumber[index] not in tailnumber_dict:
        tailnumber_dict[equip] = [] 
        print('Tailnumber does not exist for: ', equip, tailnumber[index])
    else:
        print('Tailnumber already exists for: ', equip, tailnumber[index])

    if equip not in equip_count_dict:
        equip_count_dict[equip] = 0

    if equip_count_dict[equip] >= 5:
        print('Already 5 inversions for: ', equip, equip_count_dict[equip])

    else:
        print('This ' + str(equip), ' has ' + str(equip_count_dict[equip]) + ' inversions')
    
    for i in range(len(stations)):
        if stations[i] == sta:
            seismo_lat = seismo_latitudes[i]
            seismo_lon = seismo_longitudes[i]
            elev = elevations[i]
            break
    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)

    if rerun_fig == False:
        output = open('output/' + equip + 'data_atmosphere_full.csv', 'a')

    input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'
    
    try:
        file =  open(input_files, 'r') 
    except:
        #print('No file for: ', date, flight_num, sta)

        output = str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'
        ht_c = datetime.fromtimestamp(closest_time, tz=timezone.utc)
        h_c = ht_c.hour
        comand = f'ncpag2s.py point --date {date[0:4]}-{date[4:6]}-{date[6:8]} --hour {h_c} --lat {lat} --lon {lon} --output {output}'

        os.system(comand)
        try:
            file =  open(output, 'r') 
        except:
            print(output)
