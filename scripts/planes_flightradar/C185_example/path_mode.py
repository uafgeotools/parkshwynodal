import pyproj
import pygmt
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Ensure repository root is on sys.path so local package 'src' can be imported
repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
from src.doppler_funcs import find_closest_point

folder_path = '/home/irseppi/REPOSITORIES/'
#Load seismometer data
seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/parkshwy_nodes.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']

# Define the UTM projection with zone 6 and WGS84 ellipsoid
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

#Convert to UTM coordinates
seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]


flight_num = 529754214
tail_num = 10512184

file = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/output/NGT_flight_param_inv_DB.txt', sep=",")
flight_nums = file['Flight_Number']
freq_peaks = file['Meas_Source_Frequency_Array']
nodes = file['Station']
x_airport, y_airport = utm_proj(-150.1072713049972,62.30091781635389)

# Create a dictionary to store the color for each tail number
all_med = {}
points_lat = {}
points_lon = {}
flight_lat = {}
flight_lon = {}
flight_alt = {}
UTM_km_x = {}
UTM_km_y = {}
points_sta_lat = {}
points_sta_lon = {}


all_med[flight_num] = []
points_lat[flight_num] = []
points_lon[flight_num] = []
flight_lat[flight_num] = []
flight_lon[flight_num] = []
flight_alt[flight_num] = []
UTM_km_x[flight_num] = []
UTM_km_y[flight_num] = []
points_sta_lat[flight_num] = []
points_sta_lon[flight_num] = []


flight_latitudes = []
flight_longitudes = []
alt = []

flight_file = '/scratch/irseppi/nodal_data/flightradar24/20190221_positions/20190221_529754214.csv'
flight_data = pd.read_csv(flight_file, sep=",")
flight_latitudes.extend(list(flight_data['latitude']))
flight_longitudes.extend(list(flight_data['longitude']))
alt_ft = list(flight_data['altitude'])
alt.extend([a * 0.3048 for a in alt_ft])

# Convert flight latitude and longitude to UTM coordinates
flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
flight_utm_x, flight_utm_y = zip(*flight_utm)

# Convert UTM coordinates to kilometers
flight_utm_x_km = [x / 1000 for x in flight_utm_x]
flight_utm_y_km = [y / 1000 for y in flight_utm_y]
flight_path = [(x,y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]

flight_lat[flight_num].extend(flight_latitudes)
flight_lon[flight_num].extend(flight_longitudes)
UTM_km_x[flight_num].extend(flight_utm_x_km)
UTM_km_y[flight_num].extend(flight_utm_y_km)
flight_alt[flight_num].extend(alt)

for i,f_num in enumerate(flight_nums):
    if int(f_num) != int(flight_num):
        continue
    node = int(nodes[i])
    # Find find the data for the example station
    for s in range(len(seismo_data)):
        if str(node) == str(stations[s]):
            seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s]) 
            closest_station_lat = seismo_latitudes[s]
            closest_station_lon = seismo_longitudes[s]
            break
        else:
            continue
    
    peaks = np.array(freq_peaks[i].strip('[]').split(' '), dtype=float)

    closest_p, dist_km, index = find_closest_point(flight_path, seismometer)
    closest_x, closest_y = closest_p
    closest_lon, closest_lat = utm_proj(closest_x * 1000, closest_y * 1000, inverse=True)

    f1 = []
    peak_old = 0
    for peak in peaks:
        if np.abs(float(peak) - float(peak_old)) < 10:
            continue
        if peak == peaks[0]:
            peak_old = peak
            continue

        diff = float(peak) - float(peak_old)
        if diff > 22 or diff < 18:
            peak_old = peak
            continue
        f1.append(diff)
        peak_old = float(peak)
    if len(f1) < 1:
        continue

    all_med[flight_num].extend([np.nanmedian(f1)])
    points_lat[flight_num].extend([closest_lat])
    points_lon[flight_num].extend([closest_lon])
    points_sta_lat[flight_num].extend([closest_station_lat])
    points_sta_lon[flight_num].extend([closest_station_lon])


f_lat = flight_lat[flight_num]
f_lon = flight_lon[flight_num]
alt_t = flight_alt[flight_num]
lat = np.array(points_lat[flight_num])
lon = np.array(points_lon[flight_num])
med = np.array(all_med[flight_num])
fxx = np.array(UTM_km_x[flight_num])
fyy = np.array(UTM_km_y[flight_num])
p_sta_lat = np.array(points_sta_lat[flight_num])
p_sta_lon = np.array(points_sta_lon[flight_num])
P_lat = np.array(points_lat[flight_num])
P_lon = np.array(points_lon[flight_num])

P_utm = [utm_proj(lon, lat) for lat, lon in zip(P_lat, P_lon)]
P_utm_x, P_utm_y = zip(*P_utm)


P_utm_x_km = [x / 1000 for x in P_utm_x]
P_utm_y_km = [y / 1000 for y in P_utm_y]
dist_point = []
elev_point = []

#Create array of distance along path 
dist_p = np.zeros(len(fxx))
dist_hold = 0
for i in range(1,len(fxx)): 
    dist_p[i] = np.sqrt((np.array(fxx)[i] - np.array(fxx)[i-1]) ** 2 + (np.array(fyy)[i] - np.array(fyy)[i-1]) ** 2) + dist_hold
    dist_hold = dist_p[i]


for j in range(len(P_utm_x_km)):
    for i in range(1,len(fxx)+1):
        # Check if the point (P_utm_x_km[j], P_utm_y_km[j]) lies between the two points (fxx[i], fyy[i]) and (fxx[i+1], fyy[i+1])
        if min(np.array(fxx)[i-1], np.array(fxx)[i]) <= P_utm_x_km[j] <= max(np.array(fxx)[i-1], np.array(fxx)[i]) and \
        min(np.array(fyy)[i-1], np.array(fyy)[i]) <= P_utm_y_km[j] <= max(np.array(fyy)[i-1], np.array(fyy)[i]):
            # Calculate the distance along the path to the point (P_utm_x_km[j], P_utm_y_km[j])
            segment_length = np.sqrt((np.array(fxx)[i] - np.array(fxx)[i-1]) ** 2 + (np.array(fyy)[i] - np.array(fyy)[i-1]) ** 2)
            projection_factor = np.sqrt((P_utm_x_km[j] - np.array(fxx)[i-1]) ** 2 + (P_utm_y_km[j] - np.array(fyy)[i-1]) ** 2) / segment_length
            projected_dist = dist_p[i-1] + projection_factor * segment_length
            dist_point.append(projected_dist)
                
            # Interpolate the altitude at the point (P_utm_x_km[j], P_utm_y_km[j])
            interpolated_alt = alt_t[i-1] + projection_factor * (alt_t[i] - alt_t[i-1])
            elev_point.append(interpolated_alt)
            break

        else:
            continue

fig = pygmt.Figure()
with pygmt.config(MAP_DEGREE_SYMBOL= "none"):
    with fig.subplot(
        nrows=1,
        ncols=2,
        figsize=("16c", "25c"),
        margins=["0.1c", "0.1c"],autolabel=False,
    ): 
        with fig.set_panel(panel=[0,1]):
            grid = pygmt.datasets.load_earth_relief(resolution="15s", region=[-151.2, -150.05, 62.29, 63.15], registration="pixel")
            with pygmt.config(MAP_FRAME_TYPE='plain', FORMAT_GEO_MAP="ddd.xx"):
                proj = "M15c"
                fig.grdimage(grid=grid, projection=proj, frame="a", cmap="geo")
                fig.colorbar(frame=["a1000", "x+lElevation, m"], position="JMR+o8c/6c+w11.5c/0.5c")
                fig.plot(x=np.array(f_lon), y=np.array(f_lat), pen="1p,black", projection=proj)
                for i in range(len(f_lat) - 1):
                    if i == 4:
                        angle = np.arctan2(np.array(f_lat)[i + 4] - np.array(f_lat)[i], np.array(f_lon)[i + 4] - np.array(f_lon)[i])
                        angle = np.degrees(angle)
                        fig.plot(x=[np.array(f_lon)[i]], y=[np.array(f_lat)[i]], style="v0.7c+e", direction=[[angle-22], [0.7]], fill='black', pen="1p,black", region=[-151.2, -150.05, 62.29, 63.15],projection=proj)
                    elif i == len(f_lat) - 8:
                        angle = np.arctan2(np.array(f_lat)[i + 2] - np.array(f_lat)[i], np.array(f_lon)[i + 2] - np.array(f_lon)[i])
                        angle = np.degrees(angle)

                        fig.plot(x=[np.array(f_lon)[i]], y=[np.array(f_lat)[i]], style="v0.7c+e", direction=[[angle-21], [1.5]], fill='black', pen="1p,black", region=[-151.2, -150.05, 62.29, 63.15],projection=proj)
                    else:
                        continue
                fig.plot(x=seismo_longitudes, y=seismo_latitudes, style="x0.2c", pen="01p,black", projection=proj)
                fig.plot(x=-150.1072713049972, y=62.30091781635389, style="x0.3c", pen="02p,pink", projection=proj)
                pygmt.makecpt(cmap="gmt/seis", series=[np.min(med)-0.01, np.max(med)+0.01])
                yy = fig.plot(x=lon, y=lat, style="c0.3c", fill=med, pen="black", cmap=True, projection=proj)
                fig.colorbar(frame=["a0.5f0.1", 'xaf+l\u0394'+'f, Hz'], position="JMR+o8c/-6.5c+w11.5c/0.5c")

                zoom_region = [np.min(lon) - 0.01, np.max(lon) + 0.01, np.min(lat) - 0.01, np.max(lat) + 0.01]
                rectangle = [[zoom_region[0], zoom_region[2], zoom_region[1], zoom_region[3]]]
                fig.plot(data=rectangle, style="r+s", pen="0.5p,black", projection=proj)
                # add a km bar scale
                fig.basemap(region=[-151.2, -150.05, 62.29, 63.15], projection=proj, map_scale="jBL+w5k+f+o0.5c/0.5c")

                fig.text(
                    text="b)",
                    x=-151.19,
                    y=63.13,
                    font="14p,Helvetica-Bold,black",
                    fill="white",
                    justify="LB",
                    projection=proj,
                )
        with fig.set_panel(panel=[0,0]):
            zoom_region = [np.min(lon) - 0.01, np.max(lon) + 0.01, np.min(lat) - 0.005, np.max(lat)+ 0.015]
            with pygmt.config(MAP_FRAME_TYPE = 'plain',FORMAT_GEO_MAP="ddd.xx"):
                grid_inset = pygmt.datasets.load_earth_relief(resolution="15s", region=zoom_region, registration="pixel")
                cent_lon =  str(((np.max(lon) + 0.01) - np.min(lon) - 0.01)/2 + (np.min(lon) - 0.01))
                cent_lat = str(((np.max(lat) + 0.01) - np.min(lat) - 0.01)/2 + (np.min(lat) - 0.01))
                proj = "M6.8c"
                cmap_limits = [float(np.min(grid)), float(np.max(grid))]  # Get min and max elevation values
                pygmt.makecpt(cmap="geo", series=cmap_limits, continuous=True)
            
                fig.grdimage(grid=grid_inset, region=zoom_region, projection=proj,frame="a", cmap=True)
                fig.plot(x=np.array(f_lon), y=np.array(f_lat), projection=proj, pen="1p,black") 

                for i in range(len(f_lat) - 1):
                    if i == 26:
                        angle = np.arctan2(np.array(f_lat)[i + 4] - np.array(f_lat)[i], np.array(f_lon)[i + 4] - np.array(f_lon)[i])
                        angle = np.degrees(angle)
                        fig.plot(x=[np.array(f_lon)[i]], y=[np.array(f_lat)[i]], style="v0.9c+e", direction=[[angle-21], [1]], fill='black', pen="1p,black", region=zoom_region,projection=proj)
                    elif i == len(f_lat) - 18:
                        angle = np.arctan2(np.array(f_lat)[i + 10] - np.array(f_lat)[i], np.array(f_lon)[i + 10] - np.array(f_lon)[i])
                        angle = np.degrees(angle)

                        fig.plot(x=[np.array(f_lon)[i]], y=[np.array(f_lat)[i]], style="v0.9c+e", direction=[[angle-16], [1]], fill='black', pen="1p,black", region=zoom_region,projection=proj)
                    else:
                        continue
                
                for tt in range(len(med)):
                    x = [lon[tt], p_sta_lon[tt]]
                    y = [lat[tt], p_sta_lat[tt]]
                    fig.plot(x=x, y=y, pen="1p,black,-", connection="n", projection=proj)

                fig.plot(x=seismo_longitudes, y=seismo_latitudes, projection=proj, style="x0.6c", pen="2p,black")
                pygmt.makecpt(cmap="gmt/seis", series=[np.min(med)-0.01, np.max(med)+0.01]) 
                yy = fig.plot(x=lon, y=lat, style="c0.4c", fill=med, projection=proj, pen="black", cmap=True) 
                fig.basemap(region=zoom_region, projection=proj, map_scale="jBL+w1k+f+o0.5c/0.5c")
                fig.text(
                    text="a)",
                    x=zoom_region[0] + 0.0025,
                    y=zoom_region[3] - 0.0014,
                    font="14p,Helvetica-Bold,black",
                    fill="white",
                    justify="LT",
                    projection=proj,
                )
    fig.shift_origin(yshift="-8c")
    proj = "X24.5c/6c"
    with fig.subplot(nrows=1, ncols=1, figsize=("27c", "10c"), margins=["0.1c", "0.1c"],autolabel=False):

        points = pd.DataFrame({'longitude': [], 'latitude': []})
        for i in range(1, len(f_lon)):
            lon_segment = np.linspace(f_lon[i-1], f_lon[i], 4)  # 4 intermediate points + start and end
            lat_segment = np.linspace(f_lat[i-1], f_lat[i], 4)
            points = pd.concat([points, pd.DataFrame({'longitude': lon_segment, 'latitude': lat_segment})], ignore_index=True)

        elevation_data = pygmt.grdtrack(
            grid=grid,
            points=points,
            newcolname="elevation",
            resample="r",
        )      

        ev = elevation_data['elevation'].values  
        ev_lat = elevation_data['latitude'].values
        ev_lon = elevation_data['longitude'].values

        # Convert flight latitude and longitude to UTM coordinates
        ev_utm = [utm_proj(lon, lat) for lat, lon in zip(ev_lat, ev_lon)]
        ev_utm_x, ev_utm_y = zip(*ev_utm)

        ev_utm_x_km = [x / 1000 for x in ev_utm_x]
        ev_utm_y_km = [y / 1000 for y in ev_utm_y]

        interpolated_dist_p = np.zeros(len(ev_utm_x_km))
        dist_hold = 0
        for i in range(1,len(ev_utm_x_km)):
            interpolated_dist_p[i] = np.sqrt((np.array(ev_utm_x_km)[i] - np.array(ev_utm_x_km)[i-1]) ** 2 + (np.array(ev_utm_y_km)[i] - np.array(ev_utm_y_km)[i-1]) ** 2) + dist_hold
            dist_hold = interpolated_dist_p[i]
        # Interpolate elevation for every 1 km along the path
        interp_km = np.arange(0, np.max(interpolated_dist_p), 0.1)
        interp_elev = np.interp(interp_km, interpolated_dist_p, ev)

        # Now interp_km contains every 1 km along the path, interp_elev contains the corresponding elevation
        # You can use interp_km and interp_elev for further plotting or analysis
        distance_grid, elevation_grid = np.meshgrid(
            interp_km,
            np.linspace(0, np.max(alt_t) + 100, len(interp_km))
        )

        # Fill the mesh grid with elevation values for color mapping
        color_fill = np.full_like(distance_grid, 0)  # Initialize with NaN
        for row in range(len(elevation_grid)):
            for col in range(len(elevation_grid[row])):
                if elevation_grid[row, col] <= float(interp_elev[col]):
                    color_fill[row, col] = elevation_grid[row, col]
                    if elevation_grid[row, col] <= np.min(interp_elev):
                        color_fill[row, col] = np.min(interp_elev)
                else:
                    color_fill[row, col] = np.nan  


        # Turn grids into 1D arrays
        distance_grid = distance_grid.flatten()
        elevation_grid = elevation_grid.flatten()
        color_fill = color_fill.flatten()

        num_cells = 600  # between 600 and 650
        space_x = (np.max(distance_grid) - np.min(distance_grid)) / num_cells
        space_y = (np.max(elevation_grid) - np.min(elevation_grid)) / num_cells
        grid_y_min, grid_y_max = np.min(elevation_grid), np.max(elevation_grid)
        outline_y_min, outline_y_max = np.min(ev), np.max(ev)

        # Choose the larger range to ensure both the grid and outline fit
        y_min = min(grid_y_min, outline_y_min)
        y_max = max(grid_y_max, outline_y_max)

        # Ensure the region covers both grid and outline
        prof_region = [0, np.max(distance_grid), y_min, y_max]
        fig.basemap(
            region=prof_region,  # x_min, x_max, y_min, y_max
            projection=proj,
            frame=["WSrt", "xa20+lDistance, km", "ya1000+lElevation, m"], 
        )

        fig.plot(
            x=[0, np.max(dist_p), np.max(dist_p), 0],
            y=[0, 0, np.max(alt_t)+100, np.max(alt_t)+100],
            fill="lightblue",
            projection=proj,
            close=True
        )

        c_fill = pygmt.xyz2grd(
            x=distance_grid, 
            y=elevation_grid, 
            z=color_fill, 
            projection=proj, 
            region=prof_region, 
            spacing=(space_x, space_y)
        )
        cmap_limits = [float(np.min(grid)), float(np.max(grid))]  # Include light blue value in the range
        pygmt.makecpt(cmap="geo", series=cmap_limits, continuous=True)

        # Apply the color palette table to the grid
        fig.grdimage(grid=c_fill, projection=proj, region=prof_region, cmap=True, nan_transparent=True)
        fig.plot(x=np.array(dist_p), y=np.array(alt_t), pen="1p,black", region=prof_region, projection=proj)
        fig.plot(
            x=np.array(interpolated_dist_p),
            y=np.array(ev),
            pen="2p,black",
            projection=proj,
            region=prof_region,
        )

        pygmt.makecpt(cmap="gmt/seis", series=[np.min(med)-0.01, np.max(med)+0.01]) 

        fig.plot(
            x=np.array(dist_point),
            y=np.array(elev_point),
            style="c0.3c", 
            fill=med, 
            pen="black", 
            cmap=True,
            projection=proj,
            region=prof_region,
        )

        fig.text(
            text="c)",
            x=1.5,
            y=np.max(prof_region) - 280,
            font="14p,Helvetica-Bold,black",
            fill="white",
            justify="LB",
            projection=proj,
        )

fig.savefig(folder_path + "flight_path_10512184.pdf", dpi=300)
fig.show(verbose="i")
