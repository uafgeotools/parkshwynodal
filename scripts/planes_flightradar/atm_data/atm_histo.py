import numpy as np
import pandas as pd
import os
import sys
import json
from matplotlib import pyplot as plt
from datetime import datetime, timezone
from pyproj import Proj
from pathlib import Path
# Ensure repository root is on sys.path so local package 'src' can be imported
repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
from src.doppler_funcs import speed_of_sound, add_wind_vector

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/parkshwy_nodes.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sensor_elevations = seismo_data['Elevation']
sensors = seismo_data['Station']

air_temp_array = []
air_wind_array = []
air_c_array = []
sta_temp_array = []
sta_wind_array = []
sta_c_array = []
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
for li in file_in.readlines():
    text = li.split(',')
    flight_num = text[1]
    date = text[0]
    sta = text[9]
    time = float(text[5])
    start_time = time - 120
    equip = text[10]
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
            print('No file for: ', date, flight_num, sta)
            continue
    else:
        print('No file for: ', date, flight_num, sta)
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
                elif item['parameter'] == 'V':
                    meridional_winds = float(item['values'][z_index])
                elif item['parameter'] == 'U':
                    zonal_winds = float(item['values'][z_index])

            wind_air, az_air = add_wind_vector(zonal_winds, meridional_winds)
            c_air = speed_of_sound(Tc_air)
            air_temp_array.append(Tc_air)
            air_c_array.append(c_air)
            #Figure out component of wind in direction of sensor (90 degrees to heading direction?)
            air_wind_array.append(abs(wind_air))
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
                    elif item['parameter'] == 'V':
                        meridional_winds = float(item['values'][z_index])
                    elif item['parameter'] == 'U':
                        zonal_winds = float(item['values'][z_index])
            wind_sta, az_sta = add_wind_vector(zonal_winds, meridional_winds)
            c_sta = speed_of_sound(Tc_sta)
            sta_temp_array.append(Tc_sta)
            sta_c_array.append(c_sta)
            sta_wind_array.append(abs(wind_sta))
            file.close()

fig, ax = plt.subplots(4,3, figsize=(15, 15), sharex=False, sharey=False)

ax[0, 0].hist(air_temp_array, bins=20, color='black', alpha=0.5, edgecolor='black')
median_temp = np.median(air_temp_array)
ax[0, 0].axvline(median_temp, color='r', linestyle='--', label=str(median_temp))
ax[0, 0].set_title('Temperature Distribution (\u00b0C)')
ax[0, 0].set_ylabel('Aircraft Location Data\n' + str(len(air_temp_array)) + ' samples')
ax[0, 0].xaxis.set_label_position('top')
ax[0, 0].set_xlabel('Median aircraft temperature: {:.2f}'.format(median_temp))

ax[0, 1].hist(air_c_array, bins=20, color='black', alpha=0.5, edgecolor='black')
median_c = np.median(air_c_array)
ax[0, 1].axvline(median_c, color='r', linestyle='--', label=str(median_c))
ax[0, 1].set_title('Sound Speed Distribution (m/s)')
ax[0, 1].xaxis.set_label_position('top')
ax[0, 1].set_xlabel('Median aircraft speed of sound: {:.2f}'.format(median_c))

ax[0, 2].hist(air_wind_array, bins=20, color='black', alpha=0.5, edgecolor='black')
median_wind = np.median(air_wind_array)
ax[0, 2].axvline(median_wind, color='r', linestyle='--', label=str(median_wind))
ax[0, 2].set_title('Wind Speed Distribution (m/s)')
ax[0, 2].xaxis.set_label_position('top')
ax[0, 2].set_xlabel('Median aircraft wind speed: {:.2f}'.format(median_wind))

ax[1, 0].hist(sta_temp_array, bins=20, color='black', alpha=0.5, edgecolor='black')
median_temp = np.median(sta_temp_array)
ax[1, 0].axvline(median_temp, color='r', linestyle='--', label=str(median_temp))
ax[1, 0].set_ylabel('Sensor Location Data\n' + str(len(sta_temp_array)) + ' samples')
ax[1, 0].xaxis.set_label_position('top')
ax[1, 0].set_xlabel('Median sensor temperature: {:.2f}'.format(median_temp))

ax[1, 1].hist(sta_c_array, bins=20, color='black', alpha=0.5, edgecolor='black')
median_c = np.median(sta_c_array)
ax[1, 1].axvline(median_c, color='r', linestyle='--', label=str(median_c))
ax[1, 1].xaxis.set_label_position('top')
ax[1, 1].set_xlabel('Median sensor speed of sound: {:.2f}'.format(median_c))
ax[1, 1].set_xticks(np.arange(320,  340, 5))

ax[1, 2].hist(sta_wind_array, bins=20, color='black', alpha=0.5, edgecolor='black')
median_wind = np.median(sta_wind_array)
ax[1, 2].axvline(median_wind, color='r', linestyle='--', label=str(median_wind))
ax[1, 2].xaxis.set_label_position('top')
ax[1, 2].set_xlabel('Median sensor wind speed: {:.2f}'.format(median_wind))

ax[2, 0].hist(np.array(air_temp_array)-np.array(sta_temp_array), bins=20, color='black', alpha=0.5, edgecolor='black')
median_temp = np.median(np.array(air_temp_array)-np.array(sta_temp_array))
ax[2, 0].axvline(median_temp, color='r', linestyle='--', label=str(median_temp))
ax[2, 0].set_ylabel('Difference (Aircraft - Sensor)\n' + str(len(air_temp_array)) + ' samples')
ax[2, 0].xaxis.set_label_position('top')
ax[2, 0].set_xlabel('Median temperature difference: {:.2f}'.format(median_temp))

ax[2, 1].hist(np.array(air_c_array) - np.array(sta_c_array), bins=20, color='black', alpha=0.5, edgecolor='black')
median_c = np.median(np.array(air_c_array) - np.array(sta_c_array))
ax[2, 1].axvline(median_c, color='r', linestyle='--', label=str(median_c))
ax[2, 1].xaxis.set_label_position('top')
ax[2, 1].set_xlabel('Median speed of sound difference: {:.2f}'.format(median_c))

ax[3, 0].hist((np.array(air_temp_array)+np.array(sta_temp_array))/2, bins=20, color='black', alpha=0.5, edgecolor='black')
median_temp = np.median((np.array(air_temp_array)+np.array(sta_temp_array))/2)
ax[3, 0].axvline(median_temp, color='r', linestyle='--', label=str(median_temp))
ax[3, 0].set_ylabel('Mean (Aircraft + Sensor)/2\n' + str(len(air_temp_array)) + ' samples')
ax[3, 0].xaxis.set_label_position('top')
ax[3, 0].set_xlabel('Median average temperature: {:.2f}'.format(median_temp))


ax[3, 1].hist((np.array(air_c_array) + np.array(sta_c_array))/2, bins=20, color='black', alpha=0.5, edgecolor='black')
median_c = np.median((np.array(air_c_array) + np.array(sta_c_array))/2)
ax[3, 1].axvline(median_c, color='r', linestyle='--', label=str(median_c))
ax[3, 1].xaxis.set_label_position('top')
ax[3, 1].set_xlabel('Median average speed of: {:.2f}'.format(median_c))

# Remove outline for axes with no data
for (row, col) in [(2,2), (3,2)]:
    for name, spine in ax[row, col].spines.items():
        spine.set_visible(False)
    ax[row, col].set_xticks([])
    ax[row, col].set_yticks([])
    ax[row, col].set_facecolor('none')  # Remove plot background
    plt.setp(ax[row, col].get_xticklabels(), visible=False)
    ax[row, col].tick_params(axis='x', which='both', length=0, labelbottom=False)
fig.savefig("atm_histo.pdf", dpi=500)
plt.show()
