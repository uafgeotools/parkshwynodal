from datetime import datetime, timezone
import os
from pyproj import Proj

sta_f = open('/home/irseppi/REPOSITORIES/parkshwynodal/node_crossings_db_UTM.txt','r')

# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
	text = line.split(',')
	d = text[0]
	utm_zone = '6' # Replace with the UTM zone of your coordinates
	x =  text[2]  # Replace with your UTM x-coordinate
	y = text[3]  # Replace with your UTM y-coordinate

	# Define the UTM projection
	utm_proj = Proj(proj='utm', zone=utm_zone, ellps='WGS84')

	# Convert UTM coordinates to latitude and longitude
	lon, lat = utm_proj(x, y, inverse=True)

	# Print the converted latitude and longitude
	time = float(text[5])
	ht = datetime.fromtimestamp(time, tz=timezone.utc)
	h = ht.hour

	output = str(time) + '_' + str(lat) + '_' + str(lon) + '.dat' 

	comand = f'ncpag2s.py point --date {d[0:4]}-{d[4:6]}-{d[6:8]} --hour {h} --lat {lat} --lon {lon} --output {output}'

	os.system(comand)


