import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import json
import datetime
from datetime import datetime, timezone
from pyproj import Proj
from prelude import *
from plot_func import *

def speed_of_sound(Tc):
    #Tc is the temperature in degrees celsius
    #gama = 1.4 #typical adiabatic index for air
    #c = np.sqrt(gama*R*T/M)
    c = 331.3+0.6*Tc
    return c
utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

temp_array = []
c_array = []
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
for li in file_in.readlines():
    text = li.split(',')
    flight_num = text[1]
    date = text[0]
    sta = text[9]
    time = float(text[5])
    start_time = time - 120

    # Print the converted latitude and longitude
    ht = datetime.fromtimestamp(time, tz=timezone.utc)
    h = ht.hour

    alt = float(text[4])*0.0003048 #convert between feet and km
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)

    input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'
    try:
        file =  open(input_files, 'r') #as file:
    except:
        print('No file for: ', date, flight_num, sta)
        continue
    data = json.load(file)

    # Extract metadata
    metadata = data['metadata']
    sourcefile = metadata['sourcefile']
    datetim = metadata['time']['datetime']
    latitude = metadata['location']['latitude']
    longitude = metadata['location']['longitude']
    parameters = metadata['parameters']

    # Extract data
    data_list = data['data']

    # Convert data to a DataFrame
    data_frame = pd.DataFrame(data_list)

    # Find the "Z" parameter and extract the value at index
    z_index = None
    hold = np.inf
    for item in data_list:
        if item['parameter'] == 'Z':
            for i in range(len(item['values'])):
                if abs(float(item['values'][i]) - float(alt)) < hold:
                    hold = abs(float(item['values'][i]) - float(alt))
                    z_index = i

    for item in data_list:
        if item['parameter'] == 'T':
            Tc = - 273.15 + float(item['values'][z_index])
            temp = Tc
            Tc_low = - 273.15 + float(item['values'][1])
    c = speed_of_sound(Tc)
    c_low = speed_of_sound(Tc_low)
    sound_speed = c
    print(f"Speed of sound: {c} m/s")
    print(f"Speed of sound low: {c_low} m/s")
    print(f"Temperature: {Tc} degrees Celsius")
    print(f"Temperature_low: {Tc_low} degrees Celsius")
    temp_array.append(Tc)
    c_array.append(c)
    file.close()
plt.figure()
plt.hist(temp_array, bins=20)
median_temp = np.median(temp_array)
plt.axvline(median_temp, color='r', linestyle='--', label=str(median_temp))
plt.legend()
plt.show()

plt.figure()
plt.hist(c_array,bins=20)
median_c = np.median(c_array)
plt.axvline(median_c, color='r', linestyle='--', label=str(median_c))
plt.legend()
plt.show()
