import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import json
import datetime
from datetime import datetime, timezone
from pyproj import Proj
from src.doppler_funcs import speed_of_sound

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sensor_elevations = seismo_data['Elevation']
sensors = seismo_data['Station']

est_factor = []
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
for li in file_in.readlines():
    text = li.split(',')
    flight_num = text[1]
    date = text[0]
    sta = text[9]
    time = float(text[5])
    equip = text[10]
    v = float(text[7])
    if equip in ['B737', 'B738', 'B739', 'B77W', 'B772', 'B788', 'B789', 'B763', 'B744','B733','B732','B77L','B748','CRJ2', 'A332', 'A359', 'E75S']:
        continue
    for i, sensor in enumerate(sensors):
        if str(sensor) == str(sta):
            index = i
            break

    sta_lat = seismo_latitudes[index]
    sta_lon = seismo_longitudes[index]

    # Print the converted latitude and longitude
    ht = datetime.fromtimestamp(time, tz=timezone.utc)
    h = ht.hour

    alt_m = float(text[6]) 
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)

    input_files = ['/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat','/scratch/irseppi/nodal_data/plane_info/atmosphere_data_nodes/' + str(time) + '_' + str(sta_lat) + '_' + str(sta_lon) + '.dat']

    if os.path.exists(input_files[0]):
        if os.path.exists(input_files[1]):
            pass
        else:
            continue
    else:
        continue
    for ii,input_file in enumerate(input_files):
        file = open(input_file, 'r')
        data = json.load(file)

        # Extract metadata
        metadata = data['metadata']
        parameters = metadata['parameters']

        # Extract data
        data_list = data['data']

        # Convert data to a DataFrame
        data_frame = pd.DataFrame(data_list)

        if ii == 0:
            # Find the "Z" parameter and extract the value at index
            z_index = None
            hold = np.inf
            for item in data_list:
                if item['parameter'] == 'Z':
                    for i in range(len(item['values'])):
                        if abs(float(item['values'][i]) - float(alt_m/1000)) < hold:
                            hold = abs(float(item['values'][i]) - float(alt_m/1000))
                            z_index = i

            for item in data_list:
                if item['parameter'] == 'T':
                    Tc_air = - 273.15 + float(item['values'][z_index])

            c_air = speed_of_sound(Tc_air)
            file.close()

        else:
            # Find the "T" parameter and extract the value at index
            for item in data_list:
                # Find the "Z" parameter and extract the value at index
                z_index = None
                hold = np.inf
                for item in data_list:
                    if item['parameter'] == 'Z0':
                        ground_height = float(item['values'][0])
                        break 

                for item in data_list:
                    if item['parameter'] == 'Z':
                        for i in range(len(item['values'])):
                            if abs(float(item['values'][i]) - float(ground_height)) < hold:
                                hold = abs(float(item['values'][i]) - float(ground_height))
                                z_index = i
                for item in data_list:
                    if item['parameter'] == 'T':
                        Tc_sta = - 273.15 + float(item['values'][z_index])


            c_sta = speed_of_sound(Tc_sta)
            file.close()
    c_avg = (c_air + c_sta)/2
    term = (1- (v/c_avg)**2) 
    est_factor.append(term)
plt.figure()
plt.hist(est_factor, bins=20, color='black', alpha=0.5, edgecolor='black')
median = np.median(est_factor)
mad = np.nanmedian(np.abs(est_factor - median))
plt.axvline(median, color='r', linestyle='--')
plt.axvline(median-mad, color='r', linestyle='--')
plt.axvline(median+mad, color='r', linestyle='--')
plt.ylabel('Samples: ' + str(len(est_factor)) + ' samples')
plt.title('Median f\u209B/$F_c$ Difference: {:.2f} ± {:.2f}'.format(median, mad), fontsize=10)
plt.show()
