from datetime import datetime, timezone
import os
import pandas as pd
from pyproj import Proj

seismo_data = pd.read_csv('input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']

sta_f = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
	text = line.split(',')
	d = text[0]
	sta = text[9]
	print(sta)
	#index = stations.loc[stations == sta].index
	for i, station in enumerate(stations):
		if str(station) == str(sta):
			index = i
			break

	lat = seismo_latitudes[index]
	lon = seismo_longitudes[index]

	# Print the converted latitude and longitude
	time = float(text[5])
	ht = datetime.fromtimestamp(time, tz=timezone.utc)
	h = ht.hour

	output = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data_nodes/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat' 

	comand = f'ncpag2s.py point --date {d[0:4]}-{d[4:6]}-{d[6:8]} --hour {h} --lat {lat} --lon {lon} --output {output}'
	print(comand)
	os.system(comand)


