import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.doppler_funcs import *
import os

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station_elevations = seismo_data['Elevation']
stations = seismo_data['Station']

cc_array = []
fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=False, layout='constrained')
folder_list = ['/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results/','/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results_no_g_truth/']
c_array = []
for jj,folder in enumerate(folder_list):
    for filename in os.listdir(folder):
        if filename.endswith('.txt') and filename != 'spec_error_log.txt':
            filepath = os.path.join(folder, filename)
            with open(filepath, 'r') as file:
                for line in file.readlines():
                    lines = line.strip().split(',')
                    if len(lines) > 8:
                        c_array.append(float(lines[8]))



    bin = 50 #int((np.max(c_array) - np.min(c_array)) / 3)

    axs[jj].hist(c_array, bins=bin, color='k', edgecolor='black', alpha=0.5)
    axs[jj].axvline(np.median(c_array), color='red', linestyle='--', linewidth=2)
    if jj == 0:
        axs[jj].set_title('With Ground Truth\n Median Sound Speed (m/s): {:.2f}'.format(np.median(c_array)), fontsize=14)
    else:
        axs[jj].set_title('Without Ground Truth\n Median Sound Speed (m/s): {:.2f}'.format(np.median(c_array)), fontsize=14)

air_c_array = []
sta_c_array = []
c_avg_array = []
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
for li in file_in.readlines():
    text = li.split(',')
    flight_num = text[1]
    date = text[0]
    sta = text[9]
    time = float(text[5])
    equip = text[10]
    if equip in ['B737', 'B738', 'B739', 'B77W', 'B772', 'B788', 'B789', 'B763', 'B744','B733','B732','B77L','B748','CRJ2', 'A332', 'A359', 'E75S']:
        continue
    for i, station in enumerate(stations):
        if str(station) == str(sta):
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

            c_air = speed_of_sound(Tc_air)
            air_c_array.append(c_air)
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
            sta_c_array.append(c_sta)
            file.close()
    c_avg = (c_air + c_sta) / 2
    c_avg_array.append(c_avg)     



axs[2].hist(c_avg_array, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[2].axvline(np.median(c_avg_array), color='red', linestyle='--', linewidth=2)

axs[2].set_title('G2S Model\n Median Sound Speed (m/s): {:.2f}'.format(np.median(c_avg_array)), fontsize=14)
plt.show()
