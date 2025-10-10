import numpy as np
import pandas as pd
import json
import obspy
import datetime
from datetime import datetime, timezone
from pyproj import Proj
from src.doppler_funcs import *
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

rerun_fig = True #False #Flag rerun the figures without saving the inversion results = True

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
    if equip != 'DH8A' and int(sta) != 1263:
        continue
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
        tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - 120, tr[2].stats.starttime + (mins * 60) + secs + 120)
        data = tr[2][:]
        fs = int(tr[2].stats.sampling_rate)
        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
        t_wf = tr[2].times()
        if len(data) == 0:
            data = tr[1][:]
            fs = int(tr[1].stats.sampling_rate)
            title = f'{tr[1].stats.network}.{tr[1].stats.station}.{tr[1].stats.location}.{tr[1].stats.channel} − starting {tr[1].stats["starttime"]}'                        
            t_wf = tr[1].times()
            if len(data) == 0:
                data = tr[0][:]
                fs = int(tr[0].stats.sampling_rate)
                title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
                t_wf = tr[0].times()
    except:
        try:
            p = "/scratch/irseppi/500sps/2019_0" + str(month) + "_" + str(day) + "/ZE_" + str(sta) + "_DPZ.msd"
            tr = obspy.read(p)

            tr.trim(tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs - 120, tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs + 120)

            data = tr[0][:]
            fs = int(tr[0].stats.sampling_rate)
            title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
            t_wf = tr[0].times()
        except:
            print(p)
            continue

    # Compute spectrogram
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
    # Error here with division by zero ##fix this
    spec, MDF = remove_median(Sxx)

    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    t0 = tarrive-start_time
    v0 = speed_mps
    height_m = alt - elev
    l = np.sqrt(dist_m**2 + (height_m)**2)

    tf = np.arange(0, 240, 1)
    plt.figure(figsize=(12, 6))
    plt.pcolormesh(times, frequencies, spec, cmap='pink_r',shading='gouraud', vmin=vmin, vmax=vmax)
    plt.colorbar(label='Amplitude (dB)')
    coords = doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time,make_picks=False) 
    coords_array = np.array(coords)
    try:
        fobs = coords_array[:, 1]
        tobs = coords_array[:, 0]
        plt.scatter(tobs, fobs, color='red', s=15, marker='x')
    except:
        pass
    if len(coords) == 0:
        print('No picks for: ', date, flight_num, sta)
        continue
    # Convert the list of coordinates to a numpy array
    coords_array = np.array(coords)


    peaks, freqpeak =  overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time, t0, 120, make_picks=False)
    plt.scatter(freqpeak, peaks, color='red', s=15, marker='x')

    corridor_width = 10

    peaks_assos = []
    fobs = []
    tobs = []

    for pp in range(len(peaks)):
        tprime = freqpeak[pp]
        ft0p = peaks[pp]
        f0 = calc_f0(tprime, t0, ft0p, v0, l, c)
        ft = calc_ft(times,  t0, f0, v0, l, c)
        
        maxfreq = []
        coord_inv = []
        ttt = []

        f01 = f0 + corridor_width
        f02 = f0  - corridor_width
        upper = calc_ft(times,  t0, f01, v0, l, c)
        lower = calc_ft(times,  t0, f02, v0, l, c)

        for t_f in range(len(times)):

            try:      
                tt = spec[int(np.round(lower[t_f],0)):int(np.round(upper[t_f],0)), t_f]

                #For Boeing Jets
                if equip[0] == 'B' and equip[0:1] != 'BE':
                    max_amplitude_index,_ = find_peaks(tt, prominence = 5, wlen=5, height=vmax*0.5)
                else:
                    max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
                maxa = np.argmax(tt[max_amplitude_index])
                max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa])+int(np.round(lower[t_f],0))]
                maxfreq.append(max_amplitude_frequency)
                coord_inv.append((times[t_f], max_amplitude_frequency))
                ttt.append(times[t_f])

            except:
                continue

        if len(coord_inv) > 0:
            if f0 < 200:
                coord_inv_array = np.array(coord_inv)
                mtest = [f0,v0, l, t0]
                mtest,_, F_m = invert_f(mtest, coord_inv_array, c, num_iterations=4)
                ft = calc_ft(ttt,  mtest[3], mtest[0], mtest[1], mtest[2], c)
            else:
                ft = calc_ft(ttt,  t0, f0, v0, l, c)

            delf = np.array(ft) - np.array(maxfreq)

            count = 0
            for i in range(len(delf)):
                if np.abs(delf[i]) <= (4):
                    fobs.append(maxfreq[i])
                    tobs.append(ttt[i])
                    count += 1
            peaks_assos.append(count)
    
    tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, len(peaks), peaks_assos, make_picks=True)
    plt.scatter(tobs, fobs, color='blue', s=15, marker='x')
    plt.show()
    plt.close()
