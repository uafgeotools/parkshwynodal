import pandas as pd
import pyproj

from pathlib import Path
from src.doppler_funcs import load_flight_file, closest_point_on_segment
from src.main_inv_fig_functions import plot_map

utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

# Load the seismometer location data
seismo_data = pd.read_csv('input/parkshwy_nodes.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']

# Convert latitude and longitude to UTM coordinates

seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]

seismo_utm_km = [(x, y) for x, y in zip(seismo_utm_x_km, seismo_utm_y_km)]
seismo_utm_x_km, seismo_utm_y_km = zip(*seismo_utm_km)

input = open('input/node_crossings_db_UTM.txt','r')

for line in input.readlines():
    text = line.split(',')
    date = text[0]
    month = int(date[4:6])
    day = date[6:8]
    flight_num = text[1]
    closest_x_m = float(text[2])
    closest_y_m = float(text[3])
    dist_m = float(text[4])    # Distance in meters
    dist_km = dist_m / 1000    # Convert to kilometers
    closest_time = float(text[5])
    closest_p = (closest_x_m / 1000, closest_y_m / 1000)  # Convert to kilometers
    head_avg = float(text[8])
    sta = text[9]

    file_check = '/scratch/irseppi/nodal_data/plane_info/map_all_UTM/' + str(date) + '/' + str(flight_num) + '/' + str(sta) + '/map_' + str(flight_num) + '_' + str(closest_time) + '.pdf'
    if Path(file_check).exists():
        continue

    for s, stat in enumerate(stations):
        if int(stat) == int(sta):
            seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s])
            break

    directory = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions'
    filename = date + '_' + flight_num + '.csv'
    flight_file = directory + '/' + filename
    flight_utm_x_km, flight_utm_y_km, flight_path, _, _, _, head, flight_num, date = load_flight_file(flight_file,filename)

    min_distance = float('inf')
    for i in range(len(flight_path) - 1):
        flight_utm_x1, flight_utm_y1 = flight_path[i]
        flight_utm_x2, flight_utm_y2 = flight_path[i + 1]
        x, y = seismometer
        point, d = closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, x, y)
        
        if point == None:
            continue
        elif d < min_distance:
            min_distance = d
            index = i
        else:
            continue

    print(file_check)
    

    plot_map(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km, closest_time, dist_km, index, flight_num, date, seismometer, closest_p, head, sta)


