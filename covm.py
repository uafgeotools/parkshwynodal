import numpy as np
import pandas as pd
import json
import obspy
import datetime
from datetime import datetime, timezone
from pyproj import Proj
from prelude import *
from scipy.signal import find_peaks, spectrogram
from plot_func import *
from obspy.clients.nrl import NRL
import os

nrl = NRL()

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')


# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

for li in file_in.readlines():
    text = li.split(',')
    date = text[0]

    flight_num = text[1]
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate
    dist_m = float(text[4])   # Distance in meters
    closest_time = float(text[5])
    alt = float(text[6]) 
    speed_mps = float(text[7])  # Speed in meters per second
    sta = text[9]
    equip = text[10]

    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    spec_dir = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
    
    if os.path.exists(spec_dir): # and rerun_fig == False:
        go = True
    else:
        continue
    


    for i in range(len(stations)):
        if stations[i] == sta:
            seismo_lat = seismo_latitudes[i]
            seismo_lon = seismo_longitudes[i]
            elev = elevations[i]
            break
    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)


    input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'
    
    try:
        file =  open(input_files, 'r') 
    except:
        print('No tempurature file for: ', date, flight_num, sta)
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
                if abs(float(item['values'][i]) - float(alt/1000)) < hold:
                    hold = abs(float(item['values'][i]) - float(alt/1000))
                    z_index = i

    for item in data_list:
        if item['parameter'] == 'T':
            Tc = - 273.15 + float(item['values'][z_index])

    c = speed_of_sound(Tc)
    sound_speed = c
    tarrive = calc_time(closest_time,dist_m,alt,c) 

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    time = flight_data['snapshot_id']
    timestamps = flight_data['snapshot_id']
    speed = flight_data['speed']
    altitude = flight_data['altitude']

    #Must use the tarrive time to get the correct data
    ht = datetime.fromtimestamp(tarrive, tz=timezone.utc)
    if equip[0] == 'B' and equip[0:1] != 'BE':
        wind = 120
    else:
        wind = 120
    #start_time = tarrive - wind
    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'   
    if Path(file_name).exists():
        coords = []
        if Path(file_name).is_dir():
            continue
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')

                try:
                    start_time = float(pick_data[2])
                    print('Start time: ', start_time)
                    pp = 'yep'
                    break
                except:
                    pp = 'nope'
                    break
    ht = datetime.fromtimestamp(start_time+120, tz=timezone.utc)                        
    if pp == 'nope':
        continue
    h = ht.hour
    mins = ht.minute
    secs = ht.second
    month = ht.month
    day = ht.day

    h_u = str(h+1)
    if h < 23:			
        day2 = str(day)
        if h < 10:
            h_u = '0'+str(h+1)
            h = '0'+str(h)
        else:
            h_u = str(h+1)
            h = str(h)
    else:
        h_u = '00'
        day2 = str(day+1)
    if len(str(day)) == 1:
        day = '0'+str(day)
        day2 = day

    try:
        p = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta)+".mseed"
        tr = obspy.read(p)
        tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - wind, tr[2].stats.starttime + (mins * 60) + secs + wind)
        data = tr[2][:]
        fs = int(tr[2].stats.sampling_rate)
        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
        torg = tr[2].times()
        if len(data) == 0:
            data = tr[1][:]
            fs = int(tr[1].stats.sampling_rate)
            title = f'{tr[1].stats.network}.{tr[1].stats.station}.{tr[1].stats.location}.{tr[1].stats.channel} − starting {tr[1].stats["starttime"]}'                        
            torg = tr[1].times()
            if len(data) == 0:
                data = tr[0][:]
                fs = int(tr[0].stats.sampling_rate)
                title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
                torg = tr[0].times()
    except:
        try:
            p = "/scratch/irseppi/500sps/2019_0" + str(month) + "_" + str(day) + "/ZE_" + str(sta) + "_DPZ.msd"
            tr = obspy.read(p)

            tr.trim(tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs - wind, tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs + wind)

            data = tr[0][:]
            fs = int(tr[0].stats.sampling_rate)
            title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
            torg = tr[0].times()
        except:
            print(p)
            continue
    try:
        # Compute spectrogram
        frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
        # Error here with division by zero ##fix this
        spec, MDF = remove_median(Sxx)

        middle_index =  len(times) // 2
        middle_column = spec[:, middle_index]
        vmin = 0  
        vmax = np.max(middle_column) 

        tprime0 = tarrive-start_time
        v0 = speed_mps
        height_m = alt - elev
        l = np.sqrt(dist_m**2 + (height_m)**2)

        tf = np.arange(0, 240, 1)

        coords = doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time,make_picks=False) 
        coords_array = np.array(coords)

        if len(coords) == 0:
            print('No picks for: ', date, flight_num, sta)
            continue
        # Convert the list of coordinates to a numpy array
        coords_array = np.array(coords)

        f0 = 116
        m0 = [f0, v0, l, tprime0]

        m,covm, F_m = invert_f(m0, coords_array, c, num_iterations=8)
        f0 = m[0]
        v0 = m[1]
        l = m[2]
        tprime0 = m[3]
        
        ft = calc_ft(times, tprime0, f0, v0, l, c)

        peaks = []
        p, _ = find_peaks(middle_column, distance = 7)
        corridor_width = (fs/2) / len(p) 
        if equip[0] == 'B' and equip[0:1] != 'BE':
            corridor_width = 3       
        if len(p) == 0:
            corridor_width = fs/4

        coord_inv = []

        for t_f in range(len(times)):
            upper = int(ft[t_f] + corridor_width)
            lower = int(ft[t_f] - corridor_width)
            if lower < 0:
                lower = 0
            if upper > len(frequencies):
                upper = len(frequencies)
            tt = spec[lower:upper, t_f]

            max_amplitude_index = np.argmax(tt)
            
            max_amplitude_frequency = frequencies[max_amplitude_index+lower]
            peaks.append(max_amplitude_frequency)
            coord_inv.append((times[t_f], max_amplitude_frequency))


        coord_inv_array = np.array(coord_inv)

        m,_,F_m = invert_f(m0, coord_inv_array, c, num_iterations=12)
        f0 = m[0]
        v0 = m[1]
        l = m[2]
        tprime0 = m[3]

        ft = calc_ft(times, tprime0, f0, v0, l, c)
        
        delf = np.array(ft) - np.array(peaks)
        
        new_coord_inv_array = []
        for i in range(len(delf)):
            if np.abs(delf[i]) <= 3:
                new_coord_inv_array.append(coord_inv_array[i])
        coord_inv_array = np.array(new_coord_inv_array)

        m,covm,F_m = invert_f(m0, coord_inv_array, c, num_iterations=12, sigma=5)
        #def corrcov(C):
        nx,ny = covm.shape
        if nx != ny:
            continue
            
        # c = np.sqrt(np.diag(C)).reshape(nparm,1)
        # Crho = C/(c@c.T)
        sigma = np.sqrt(np.diag(covm))
        outer_v = np.outer(sigma,sigma)
        Crho = covm / outer_v
        
        Crho[covm == 0] = 0
        #return Crho
        
        gridlines=False
        colormap='seismic'
        plt.figure(figsize=(10, 10))
        plt.imshow(Crho,cmap=colormap)
        plt.xticks(ticks=range(np.shape(Crho)[1]),labels=[str(val) for val in range(1,np.shape(Crho)[1]+1)])
        plt.yticks(ticks=range(np.shape(Crho)[0]),labels=[str(val) for val in range(1,np.shape(Crho)[0]+1)])
        if gridlines:
            xgrid = np.array(range(np.shape(Crho)[1] + 1)) - 0.5
            ygrid = np.array(range(np.shape(Crho)[0] + 1)) - 0.5
            for gridline in xgrid:
                plt.axvline(x=gridline,color='k',linewidth=1)
            for gridline in ygrid:
                plt.axhline(y=gridline,color='k',linewidth=1)
        plt.colorbar()
        plt.show()
        plt.close()
    except Exception as e:
        print(f"Error processing flight {flight_num} on {date} at station {sta}: {e}")
        continue


