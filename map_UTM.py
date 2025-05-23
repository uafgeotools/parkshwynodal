import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import pyproj
import concurrent.futures
import os
from functools import partial
from datetime import datetime, timezone
from matplotlib.patches import Rectangle
from prelude import load_flights, dist_less, make_base_dir

num_workers = os.cpu_count()

def dist_less(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km):

    f = False
    for s in range(len(flight_utm_x_km)):
        for l in range(len(seismo_utm_x_km)):
            dist_km = np.sqrt((seismo_utm_y_km[l]-flight_utm_y_km[s])**2 +(seismo_utm_x_km[l]-flight_utm_x_km[s])**2)
            if dist_km <= 2:
                f = True
                break
            else:
                continue
    return f

def closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y):
    closest_point = None
    dist_lim = np.inf

    x = [flight_utm_x1, flight_utm_x2]
    y = [flight_utm_y1, flight_utm_y2]

    if (x[1] - x[0]) == 0:
        if (y[1]-y[0]) <= 0:
            ggg = -0.001
        else:
            ggg = 0.001
        for point in np.arange(y[0], y[1], ggg):
            xx = x[0]
            yy = point
            dist_km = np.sqrt((seismo_utm_y-yy)**2 +(seismo_utm_x-xx)**2)
            
            if dist_km < dist_lim:
                dist_lim = dist_km
                closest_point = (xx,yy)
            else:
                continue

    else: 
        m = (y[1]-y[0])/(x[1]-x[0])
        b = y[0] - m*x[0]

        if (x[1] - x[0]) <= 0:
            ggg = -0.001
        else:
            ggg = 0.001
        for point in np.arange(x[0], x[1], ggg):
            xx = point
        
            yy = m*xx + b
            dist_km = np.sqrt((seismo_utm_y-yy)**2 +(seismo_utm_x-xx)**2)
            
            if dist_km < dist_lim:
                dist_lim = dist_km
                closest_point = (xx,yy)
            else:
                continue

    return closest_point, dist_lim


def find_closest_point(flight_utm, seismo_utm):
    min_distance = np.inf
    closest_point = None

    for i in range(len(flight_utm) - 1):
        flight_utm_x1, flight_utm_y1 = flight_utm[i]
        flight_utm_x2, flight_utm_y2 = flight_utm[i + 1]
        seismo_utm_x, seismo_utm_y = seismo_utm
        point, d = closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y)
        
        if point == None:
            continue
        elif d < min_distance:
            min_distance = d
            closest_point = point
            index = i
        else:
            continue
   
    return closest_point, min_distance, index
def closest_time_calc(closest_p, flight_path, timestamp, index):
    closest_x, closest_y = closest_p
    #Calculate the time of the closest point
    flight_utm_x1, flight_utm_y1 = flight_path[index]
    flight_utm_x2, flight_utm_y2 = flight_path[index + 1]

    x_timestamp_dif_vec = flight_utm_x2 - flight_utm_x1
    y_timestamp_dif_vec = flight_utm_y2 - flight_utm_y1

    cx_timestamp_dif_vec =  closest_x - flight_utm_x1
    cy_timestamp_dif_vec = closest_y - flight_utm_y1

    line_vector = (x_timestamp_dif_vec, y_timestamp_dif_vec)
    cline_vector = (cx_timestamp_dif_vec, cy_timestamp_dif_vec)

    line_magnitude = np.sqrt(line_vector[0] ** 2 + line_vector[1] ** 2)
    cline_magnitude = np.sqrt(cline_vector[0] ** 2 + cline_vector[1] ** 2)

    length_ratio = cline_magnitude / line_magnitude
    closest_time = timestamp[index] + length_ratio*(timestamp[index+1] - timestamp[index])

    return closest_time

def flight_lat_lon_to_utm(flight_latitudes, flight_longitudes):
    # Convert flight latitude and longitude to UTM coordinates
    flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
    flight_utm_x, flight_utm_y = zip(*flight_utm)

    # Convert UTM coordinates to kilometers
    flight_utm_x_km = [x / 1000 for x in flight_utm_x]
    flight_utm_y_km = [y / 1000 for y in flight_utm_y]
    flight_path = [(x,y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]
        
    return flight_utm_x_km, flight_utm_y_km, flight_path

def plot_map(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km, closest_time, d, index, flight_num, date, seismometer, closest_p, head, station):
    
    closest_x, closest_y = closest_p
    min_lon = -150.7
    max_lon = -147.3
    min_lat = 62.2
    max_lat = 65.3
    lxmin,lymin = utm_proj(min_lon, min_lat)
    lxmax, lymax = utm_proj(max_lon, max_lat)
    y = [closest_y, seismometer[1]]
    x = [closest_x, seismometer[0]]
    yy = sum(y) / len(y)
    xx = sum(x) / len(x)
    
    min_x = int(xx - 2)
    max_x = int(xx + 2)
    min_y = int(yy - 2)
    max_y = int(yy + 2)
    if d < 0.5:
        min_x = (xx - 1)
        max_x = (xx + 1)
        min_y = (yy - 1)
        max_y = (yy + 1)
    if d < 0.1:
        min_x = (xx - 0.1)
        max_x = (xx + 0.1)
        min_y = (yy - 0.1)
        max_y = (yy + 0.1)
    # Create a figure with two subplots side by side
    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 2]}) 
    fig.subplots_adjust(wspace=1)  # Adjust the spacing between subplots


    axs[0].set_xticks(np.arange(int(lxmin/1000)-7, int(lxmax/1000), 50))
    axs[0].set_yticks(np.arange(int(lymin/1000)-1, int(lymax/1000), 50))
    axs[0].set_xlabel('UTM Easting (km)')
    axs[0].set_ylabel('UTM Northing (km)')
    axs[0].set_aspect('equal')
    axs[1].set_aspect('equal')

    axs[0].grid(True, linestyle='dotted', color='gray')
    axs[1].grid(True, linestyle='dotted', color='gray')

    axs[0].scatter(seismo_utm_x_km, seismo_utm_y_km, c='#e41a1c', s = 3, label='seismometers')
    axs[0].plot(flight_utm_x_km, flight_utm_y_km, '-', c='#377eb8', lw=1, ms = 1, label='flight path')

    for i in range(1, len(flight_utm_y_km)-1, int(len(flight_utm_y_km)/5)):
        direction = np.arctan2(flight_utm_y_km[i+1] - flight_utm_y_km[i], flight_utm_x_km[i+1] - flight_utm_x_km[i])
        if (flight_utm_x_km[i+1] - flight_utm_x_km[i]) == 0:
            continue
        m = (flight_utm_y_km[i+1] - flight_utm_y_km[i])/(flight_utm_x_km[i+1] - flight_utm_x_km[i])
        if m == 0:
            continue
        b = flight_utm_y_km[i] - m*flight_utm_x_km[i]
        axs[0].quiver((flight_utm_y_km[i]-b)/m, flight_utm_y_km[i], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', headwidth = 10)

    # Set labels and title
    axs[0].set_xlim(int(lxmin/1000), int(lxmax/1000))
    axs[0].set_ylim(int(lymin/1000), int(lymax/1000))
    axs[0].tick_params(axis='both', which='major', labelsize=9)

    head_avg = (head[index]+head[index+1])/2
    converted_angle = (90 - head_avg) % 360
    heading = np.deg2rad(converted_angle)

    # Define the UTM and latitude/longitude coordinate systems
    rect = Rectangle((min_x, min_y), (max_x-min_x), (max_y-min_y), ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
    axs[0].add_patch(rect)

    axs[1].set_xticks(np.arange(min_x, max_x, np.round(((max_x - min_x) / 4), 1)))
    axs[1].set_yticks(np.arange(min_y, max_y, np.round(((max_y - min_y) / 4), 1)))

    # Draw the zoomed in map on the second subplot
    axs[1].plot(x,y, '--', c='#ff7f00')
    axs[1].plot(flight_utm_x_km, flight_utm_y_km, c='#377eb8',linestyle ='dotted')
    axs[1].scatter(flight_utm_x_km, flight_utm_y_km, c='#377eb8', s=20)

    axs[1].set_xlim(min_x, max_x)
    axs[1].set_ylim(min_y, max_y)

    axs[1].tick_params(axis='both', which='major', labelsize=9)
    axs[1].ticklabel_format(useOffset=False, style='plain')

    direction = np.arctan2(flight_utm_y_km[index+1] - flight_utm_y_km[index], flight_utm_x_km[index+1] - flight_utm_x_km[index])
        
    axs[1].quiver(closest_x, closest_y, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=8)

    axs[1].quiver(closest_x, closest_y, np.cos(heading), np.sin(heading), angles='xy', scale = 8, color='#999999')
    axs[1].scatter(closest_x, closest_y, c='#377eb8', s=50, zorder=3)
    axs[1].scatter(seismometer[0], seismometer[1], c='#e41a1c', s=50, zorder=3)

    axs[1].text(xx,yy, str(round(d, 2))+' km', fontsize=12, fontweight='bold')

    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
    con = mpatch.ConnectionPatch(xyA=(max_x, min_y), xyB=(min_x, min_y), coordsA="data", coordsB="data", axesA=axs[0], axesB=axs[1], color="black", linestyle="--")
    fig.add_artist(con)
    con = mpatch.ConnectionPatch(xyA=(max_x, max_y), xyB=(min_x, max_y), coordsA="data", coordsB="data", axesA=axs[0], axesB=axs[1], color="black", linestyle="--")
    fig.add_artist(con)
    plt.tight_layout(pad=0.2, rect=[0, 0, 1, 1])

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/map_all_UTM/' + str(date) + '/'+ str(flight_num) + '/' + str(station) + '/'
    make_base_dir(BASE_DIR)
    plt.savefig('/scratch/irseppi/nodal_data/plane_info/map_all_UTM/' + str(date) + '/' + str(flight_num) + '/' + str(station) + '/map_' + str(flight_num) + '_' + str(closest_time) + '.png', bbox_inches='tight')
    plt.close()


def load_flight_files(flight_file, filename):
    # Load flight data and extract relevant columns
    flight_data = pd.read_csv(flight_file, sep=",") 
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    timestamp = flight_data['snapshot_id']
    alt = flight_data['altitude']
    speed = flight_data['speed']
    head = flight_data['heading']
    fname = filename	
    flight_num = fname[9:18]
    date = fname[0:8]

    flight_utm_x_km, flight_utm_y_km, flight_path = flight_lat_lon_to_utm(flight_latitudes, flight_longitudes)

    return flight_utm_x_km, flight_utm_y_km, flight_path, timestamp, alt, speed, head, flight_num, date

def load_equip(date, flight_num):
    equip_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_flights.csv'
    equip_data = pd.read_csv(equip_file, sep=",")
    equip_list = equip_data['equip']
    flight_list = equip_data['flight_id']

    for i_e in range(len(equip_list)):
        if str(flight_num) == str(flight_list[i_e]) :
            equip = equip_list[i_e]
            break
    return equip

def time_check(closest_time, start_time, end_time,s_index):
    # Convert the closest time to a datetime object
    v = False
    time_airplane = datetime.fromtimestamp(closest_time, tz=timezone.utc)
    start_time_obj = datetime.strptime(start_time[s_index], "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)
    end_time_obj = datetime.strptime(end_time[s_index], "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)

    if time_airplane < start_time_obj or time_airplane > end_time_obj:
        v = True
    return v

def avg_return(alt, speed, head, index):

    alt_avg = (alt[index]+alt[index+1])/2
    alt_avg_m = alt_avg * 0.3048 #convert from feet to meters
    head_avg = (head[index]+head[index+1])/2
    speed_avg = (speed[index]+speed[index+1])/2
    speed_avg_mps = speed_avg * 0.514444 #convert from knots to meters/sec
    return alt_avg_m, speed_avg_mps, head_avg

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
    flight_utm_x_km, flight_utm_y_km, flight_path, timestamp, alt, speed, head, flight_num, date = load_flight_files(flight_file,filename)
    equip = load_equip(date, flight_num)

    # Iterate over seismometer data
    for s in range(len(seismo_data)):
        seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s])  
        station = sta[s]

        if dist_less(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km) == False:
            continue

        closest_p, d, index= find_closest_point(flight_path, seismometer)
        
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

            #output = open('input/all_station_crossing_db_UTM_updated.txt','a')
            #output.write(str(tracer))
            #output.write(str(date)+ ',' + str(flight_num) + ',' + str(closest_x_m) + ',' + str(closest_y_m) + ',' + str(dist_m) + ',' +str(closest_time) + ',' + str(alt_avg_m) + ',' + str(speed_avg_mps) + ',' + str(head_avg) + ',' + str(station) + ',' + str(equip) + ',\n')
            #output.close()
            hold_lines.append(str(date)+ ',' + str(flight_num) + ',' + str(closest_x_m) + ',' + str(closest_y_m) + ',' + str(dist_m) + ',' +str(closest_time) + ',' + str(alt_avg_m) + ',' + str(speed_avg_mps) + ',' + str(head_avg) + ',' + str(station) + ',' + str(equip) + ',\n')
            print(hold_lines)
    return hold_lines
# Load flight files and filenames
#flight_files,filenames = load_flights(2, 4, 11, 27)
flight_files,filenames = load_flights(2, 3, 10, 11)
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

tracer = [i for i in range(len(flight_files))]
with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
    output_list = executor.map(read_write_closest,flight_files,filenames,tracer)
print(output_list)
output_list = list(output_list)
output = open('input/all_station_crossing_db_UTM_updated.txt','w')
for i in range(len(output_list)):
    for j in range(len(output_list[i])):
        output.write(output_list[i][j])
output.close()
