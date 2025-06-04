import numpy as np
import pandas as pd
import pygmt
import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    radius = 6371  # Radius of earth in kilometers

    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = (math.sin(dlat / 2)**2 + 
         math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = radius * c
    return distance

colors = []
with open('input/colors.txt','r') as c_in:
	for i, line in enumerate(c_in):
		if (i + 1) % 9 == 0:
			c = str(line[0:7])
			colors.append(c)

#Load seismometer data
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
colors_dict['Unknown'] = 'magenta' 


file = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', 'r')
fig = pygmt.Figure()
with pygmt.config(MAP_DEGREE_SYMBOL= "none"):
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=[-150.5, -148.5, 62.25, 64.6], registration="pixel")
    with pygmt.config(MAP_FRAME_TYPE='plain', FORMAT_GEO_MAP="ddd.xx"):
        proj = "M15c"
        fig.grdimage(grid=grid, projection=proj, frame="a", cmap="gmt/grey")
        fig.colorbar(frame=["a1000", "x+lElevation, m"], position="JMR+o0.5c/0c+w37c/0.5c")

        flight_num_hold = 0
        x = 0
        # Iterate over each line in the file
        for lines in file.readlines():
            line = lines.split(',')
            flight_num = int(line[1])
            eq = str(line[-2])
            if eq == np.nan or eq == 'nan':
                eq = 'Unknown'

            if flight_num == flight_num_hold:
                continue
            flight_num_hold = flight_num
            flight_file = '/scratch/irseppi/nodal_data/flightradar24/'+str(line[0]) + '_positions/' + str(line[0]) + '_' + str(flight_num) + '.csv'
            flight_data = pd.read_csv(flight_file, sep=",")
            flight_lat = flight_data['latitude'] 
            flight_lon = flight_data['longitude']
            f_lat = []
            f_lon = []
            index_old = np.nan
            for i in range(len(flight_lat)):

                for t in range(len(seismo_lat)):
                    dist_km = haversine(flight_lat[i], flight_lon[i], seismo_lat[t], seismo_lon[t])
                    if dist_km < 2:
                        if abs(i-index_old) > 1:
                            if len(f_lat) > 0:
                                fig.plot(x=np.array(f_lon), y=np.array(f_lat), pen=f"0.1p,{colors_dict[eq]}", projection=proj)
                                f_lat = []
                                f_lon = []
                            else:
                                fig.plot(x=np.array(f_lon), y=np.array(f_lat), style="c0.1c", pen=f"0.1p,{colors_dict[eq]}", projection=proj)
                                f_lat = []
                                f_lon = []
                        f_lat.append(flight_lat[i])
                        f_lon.append(flight_lon[i])
                        

            if len(f_lat) > 0:
                fig.plot(x=np.array(f_lon), y=np.array(f_lat), color=colors_dict[eq], projection=proj)
            else:
                fig.plot(x=np.array(f_lon), y=np.array(f_lat), style="c0.1c", color=colors_dict[eq], projection=proj)

        fig.plot(x=seismo_lon, y=seismo_lat, style="x0.1c", pen="0.5p,black", projection=proj)

        fig.show(verbose="i")
file.close()