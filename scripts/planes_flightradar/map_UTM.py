import os
import sys
import pyproj
import concurrent.futures
import pandas as pd
from pathlib import Path

# Ensure repository root is on sys.path so local package 'src' can be imported
repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.doppler_funcs import flight_list, dist_less, time_check, find_closest_point, get_equip, closest_time_calc, load_flight_file, avg_return
from src.main_inv_fig_functions import plot_map

num_workers = os.cpu_count()

def read_write_closest(flight_file,filename,tracer):
    print((tracer/len(flight_files))*100, '%')	# Print progress
    hold_lines = []
    # Load the seismometer location data
    seismo_data = pd.read_csv('input/parkshwy_nodes.txt', sep="|")
    seismo_latitudes = seismo_data['Latitude']
    seismo_longitudes = seismo_data['Longitude']
    sta = seismo_data['Station']
    start_time = seismo_data['StartTime']
    end_time = seismo_data['EndTime']

    # Convert latitude and longitude to UTM coordinates

    seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
    seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

    # Convert UTM coordinates to kilometers
    seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
    seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]

    seismo_utm_km = [(x, y) for x, y in zip(seismo_utm_x_km, seismo_utm_y_km)]
    seismo_utm_x_km, seismo_utm_y_km = zip(*seismo_utm_km)
    flight_utm_x_km, flight_utm_y_km, flight_path, timestamp, alt, speed, head, flight_num, date = load_flight_file(flight_file,filename)
    equip = get_equip(date, flight_num)

    # Iterate over seismometer data
    for s in range(len(seismo_data)):
        seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s])  
        station = sta[s]

        if dist_less(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km) == False:
            continue

        closest_p, d, index= find_closest_point(flight_path, seismometer)
        if index == None:
            continue
        if d <= 2:

            closest_x, closest_y = closest_p
            closest_time = closest_time_calc(closest_p, flight_path, timestamp, index)

            if time_check(closest_time, start_time, end_time,s) == True:
                continue

            plot_map(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km, closest_time, d, index, flight_num, date, seismometer, closest_p, head, station)

            alt_avg_m, speed_avg_mps, head_avg = avg_return(alt, speed, head, index)

            # Convert distance to meters
            dist_m = d * 1000
            closest_x_m = closest_x * 1000
            closest_y_m = closest_y * 1000

            hold_lines.append(str(date)+ ',' + str(flight_num) + ',' + str(closest_x_m) + ',' + str(closest_y_m) + ',' + str(dist_m) + ',' +str(closest_time) + ',' + str(alt_avg_m) + ',' + str(speed_avg_mps) + ',' + str(head_avg) + ',' + str(station) + ',' + str(equip) + ',\n')
    return hold_lines

# Load flight files and filenames
flight_files,filenames = flight_list(2, 4, 11, 27)

utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

tracer = [i for i in range(len(flight_files))]
with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
    output_list = executor.map(read_write_closest,flight_files,filenames,tracer)

output_list = list(output_list)

output = open('input/all_station_crossing_db_UTM.txt','w')
for i in range(len(output_list)):
    for j in range(len(output_list[i])):
        output.write(output_list[i][j])
output.close()
