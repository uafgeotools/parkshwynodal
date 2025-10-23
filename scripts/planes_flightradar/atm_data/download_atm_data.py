from datetime import datetime, timezone
import os
import pandas as pd
from pyproj import Proj

#Created using NCPAG2S Command-Line Client: https://github.com/chetzer-ncpa/ncpag2s-clc
#Root directory of the repository must be in your PATH

output_folder_path = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/'

seismo_data = pd.read_csv('input/parkshwy_nodes.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']

# Define the UTM projection
utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')
sta_f = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

# Loop through each crossing decteted by nodes that is not a jet
for line in sta_f.readlines():
	text = line.split(',')
	d = text[0] # Date in YYYYMMDD format
	x =  text[2]  # UTM x-coordinate
	y = text[3]  # UTM y-coordinate
	sta = text[9] # Station name

	# Convert UTM coordinates to latitude and longitude
	lon, lat = utm_proj(x, y, inverse=True)

	# Extract time and hour
	time = float(text[5])
	ht = datetime.fromtimestamp(time, tz=timezone.utc)
	h = ht.hour

	output = output_folder_path + '/aircraft_loc/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat' 

	# Command to download atmospheric data at the aircraft position
	comand = f'ncpag2s.py point --date {d[0:4]}-{d[4:6]}-{d[6:8]} --hour {h} --lat {lat} --lon {lon} --output {output}'
	
	# Download atmospheric data for the given lat, lon, date, and hour
	os.system(comand)

	for i, station in enumerate(stations):
		if str(station) == str(sta):
			index = i
			break

	lat_sta = seismo_latitudes[index]
	lon_sta = seismo_longitudes[index]

	output = output_folder_path + '/nodes_loc/' + str(time) + '_' + str(lat_sta) + '_' + str(lon_sta) + '.dat' 

	# Command to download atmospheric data at the station position
	comand = f'ncpag2s.py point --date {d[0:4]}-{d[4:6]}-{d[6:8]} --hour {h} --lat {lat_sta} --lon {lon_sta} --output {output}'

	os.system(comand)


