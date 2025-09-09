import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyproj

utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

colors = []
with open('input/colors.txt','r') as c_in:
	for i, line in enumerate(c_in):
		if (i + 1) % 9 == 0:
			c = str(line[0:7])
			colors.append(c)
			
seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/nodes_stations.txt', sep="|")
seismo_lat = seismo_data['Latitude']
seismo_lon = seismo_data['Longitude']

file = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', 'r')

equip_counts = {}  # Define the "equip_counts" dictionary before the loop
for line in file:
	data = line.split(',')  # Split the line by commas
	equip = data[-2]  # Get the equipment type from the line
	if equip == np.nan or equip == 'nan':
		equip = 'Unknown'
	equip_counts[equip] = equip_counts.get(equip, 0) + 1 
file.close()

equip_counts = {k: v for k, v in sorted(equip_counts.items(), key=lambda item: item[1], reverse=True)}

colors_dict={}
for i,labels in enumerate(equip_counts.keys()):
    colors_dict[labels] = colors[i]

file = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', 'r')

flight_num_hold = 0
plt.figure(figsize=(6, 10))
# Iterate over each line in the file
for lines in file.readlines():
    line = lines.split(',')
    flight_num = int(line[1])

    eq = str(line[-2])
    x_m_UTM = float(line[2])
    y_m_UTM = float(line[3])
    
    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x_m_UTM, y_m_UTM, inverse=True)

    if eq == np.nan or eq == 'nan':
        eq = 'Unknown'

    plt.scatter(lon,lat, color = colors_dict[eq])
    plt.scatter(seismo_lon, seismo_lat, marker="x", color="black")
plt.show()
file.close()