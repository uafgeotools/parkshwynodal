import numpy as np
import pandas as pd
import numpy.linalg as la
import json
import obspy
import os
from pathlib import Path
from pyproj import Proj
from matplotlib import pyplot as plt
from datetime import datetime, timezone
from pyproj import Proj
from prelude import calc_ft, S, calc_time, speed_of_sound, calc_f0, make_base_dir, invert_f, full_inversion
from scipy.signal import find_peaks, spectrogram
from main_inv_fig_functions import doppler_picks, overtone_picks, time_picks, remove_median, plot_spectrogram, plot_spectrum
from obspy.clients.nrl import NRL


nrl = NRL()

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

rerun_fig = False #Flag rerun the figures without saving the inversion results = True

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
    
    if not os.path.exists(spec_dir): 
        continue

    for i in range(len(stations)):
        if stations[i] == sta:
            seismo_lat = seismo_latitudes[i]
            seismo_lon = seismo_longitudes[i]
            elev = elevations[i]
            break
    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)

    if rerun_fig == False:
        output = open('output/fixed_quasi/' + equip + 'data_atmosphere_full.csv', 'a')

    input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'
    
    if Path(input_files).exists():
        file =  open(input_files, 'r') 

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
    else: 
        c = 311 # Default speed of sound in m/s if no data is available
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
    
    spec_window = 120
    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'   
    if Path(file_name).exists():
        coords = []
        if Path(file_name).is_dir():
            continue
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                if len(pick_data) == 4:
                    start_time = float(pick_data[2])
                elif len(pick_data) == 3:
                    print('No start time in file: ', file_name)
                    start_time = tarrive
                    continue
                else:
                    continue
    start_time = tarrive
    if equip == 'C185':
        start_time = start_time - 120

    ht = datetime.fromtimestamp(start_time+120, tz=timezone.utc)                      

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


    waveform1 = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta)+".mseed"
    waveform2 = "/scratch/irseppi/500sps/2019_0" + str(month) + "_" + str(day) + "/ZE_" + str(sta) + "_DPZ.msd"
    if Path(waveform1).exists():
        tr = obspy.read(waveform1)
        # Trim all traces in the Stream object
        for trace in tr:
            trace.trim(trace.stats.starttime + (mins * 60) + secs - spec_window,
                   trace.stats.starttime + (mins * 60) + secs + spec_window)
        data = tr[2][:]
        fs = int(tr[2].stats.sampling_rate)
        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
        torg = tr[2].times()
        if len(data) == 0:
            data = tr[1][:]
            #tr[1].trim(tr[1].stats.starttime + (mins * 60) + secs - spec_window, tr[1].stats.starttime + (mins * 60) + secs + spec_window)
            fs = int(tr[1].stats.sampling_rate)
            title = f'{tr[1].stats.network}.{tr[1].stats.station}.{tr[1].stats.location}.{tr[1].stats.channel} − starting {tr[1].stats["starttime"]}'                        
            torg = tr[1].times()
            if len(data) == 0:
                data = tr[0][:]
                #tr[0].trim(tr[0].stats.starttime + (mins * 60) + secs - spec_window, tr[0].stats.starttime + (mins * 60) + secs + spec_window)
                fs = int(tr[0].stats.sampling_rate)
                title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
                torg = tr[0].times()
                if len(data) == 0:
                    continue
        

    elif Path(waveform2).exists():
        tr = obspy.read(waveform2)
        for trace in tr:
            trace.trim(trace.stats.starttime+ (float(h) * 3600) + (mins * 60) + secs - spec_window,
                   trace.stats.starttime + (float(h) * 3600) + (mins * 60) + secs + spec_window)
        data = tr[0][:]
        fs = int(tr[0].stats.sampling_rate)
        title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
        torg = tr[0].times()
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

    else:
        continue

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
    print(sta,equip,date,flight_num)
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    plt.show()
    plt.close()

    if len(coords) == 0:
        continue
    # Convert the list of coordinates to a numpy array
    coords_array = np.array(coords)

    f0 = (np.max(coords_array[:,1])+np.min(coords_array[:,1]))/2
    m0 = [f0, v0, l, tprime0, c]
    print('Initial guess: ', m0)

    sigma_f0 = 100
    sigma_v0 = 100
    sigma_l = 1000
    sigma_tprime0 = 200
    sigma_c = 100
    sigma_prior = [sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c]
    m,_,_, F_m = invert_f(m0,sigma_prior, coords_array, num_iterations=8)
    f0 = m[0]
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]
    
    ft = calc_ft(times, tprime0, f0, v0, l, c)

    corridor_width = 20
    if equip[0] == 'B' and equip[0:1] != 'BE':
        corridor_width = 5      

    coord_inv = []

    for t_f in range(len(times)):

        upper = int(ft[t_f] + corridor_width)
        lower = int(ft[t_f] - corridor_width)
        if lower < 0:
            lower = 0
        elif lower >= 250:
            lower = 200
        else:
            pass
        if upper > 250:
            upper = 250

        tt = spec[lower:upper, t_f]
        max_amplitude_index = np.argmax(tt)
        
        max_amplitude_frequency = frequencies[max_amplitude_index+lower]
        coord_inv.append((times[t_f], max_amplitude_frequency))

    coord_inv_array = np.array(coord_inv)

    m,_,_,F_m = invert_f(m0,sigma_prior, coord_inv_array,num_iterations=5)
    f0 = m[0]
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]

    ft = calc_ft(times, tprime0, f0, v0, l, c)
    
    delf = np.array(ft) - np.array(peaks)
    
    new_coord_inv_array = []
    for i in range(len(delf)):
        if np.abs(delf[i]) <= 3:
            new_coord_inv_array.append(coord_inv_array[i])
    coord_inv_array = np.array(new_coord_inv_array)
    m,_,_,F_m = invert_f(m0, sigma_prior, coord_inv_array, num_iterations=5, sigma=5)
        
    f0 = m[0]
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]

    peaks, freqpeak =  overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time, tprime0, tarrive, make_picks=False)

    w = len(peaks)
 
    corridor_width = (fs/2) / len(peaks) 
    if equip[0] == 'B' and equip[0:1] != 'BE':
        corridor_width = 3
    peaks_assos = []
    fobs = []
    tobs = []
    f0_array = []
    for pp in range(len(peaks)):
        tprime = freqpeak[pp]
        ft0p = peaks[pp]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        f0_array.append(f0)

        maxfreq = []
        coord_inv = []
        ttt = []

        f01 = f0 + corridor_width
        f02 = f0  - corridor_width
        upper = calc_ft(times,  tprime0, f01, v0, l, c)
        lower = calc_ft(times,  tprime0, f02, v0, l, c)

        for t_f in range(len(times)):

            if lower[t_f] < 0 or lower[t_f] > 250 or upper[t_f] > 250 or np.isnan(upper[t_f]) or np.isnan(lower[t_f]):
                continue

            tt = spec[int(np.round(lower[t_f],0)):int(np.round(upper[t_f],0)), t_f]

            #For Boeing Jets
            if str(equip[0]) == 'B' and str(equip[0:1]) != 'BE':
                max_amplitude_index,_ = find_peaks(tt, prominence = 5, wlen=5, height=vmax*0.5)
            else:
                max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
            if len(max_amplitude_index) == 0:
                continue
            maxa = np.argmax(tt[max_amplitude_index])
            max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa])+int(np.round(lower[t_f],0))]
            maxfreq.append(max_amplitude_frequency)
            coord_inv.append((times[t_f], max_amplitude_frequency))
            ttt.append(times[t_f])


        if len(coord_inv) > 0:
            if f0 < 200:
                coord_inv_array = np.array(coord_inv)
                mtest = [f0,v0, l, tprime0,c]
                mtest,_,_, F_m = invert_f(mtest,sigma_prior, coord_inv_array, num_iterations=4)
                ft = calc_ft(ttt,  mtest[3], mtest[0], mtest[1], mtest[2], mtest[4])
            else:
                ft = calc_ft(ttt,  tprime0, f0, v0, l, c)

            delf = np.array(ft) - np.array(maxfreq)

            count = 0
            for i in range(len(delf)):
                if np.abs(delf[i]) <= (4):
                    fobs.append(maxfreq[i])
                    tobs.append(ttt[i])
                    count += 1
            peaks_assos.append(count)

    tobs_hold = tobs
    if len(fobs) == 0:
        continue

    tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, w, peaks_assos, make_picks=True)

    if len(tobs) == len(tobs_hold):
        continue

    v0 = speed_mps
    height_m = alt - elev
    l = np.sqrt(dist_m**2 + (height_m)**2)    
    c = speed_of_sound(Tc)

    mprior = []
    mprior.append(v0)
    mprior.append(l)
    mprior.append(tprime0)
    mprior.append(c)
    for o in range(len(peaks_assos)):
        tprime = freqpeak[o]
        ft0p = peaks[o]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        mprior.append(float(f0))
    print("mprior:", mprior)
    plt.figure(figsize=(15, 10))
    plt.pcolormesh(times, frequencies, spec, vmin=vmin, vmax=vmax, shading='gouraud')
    plt.scatter(tobs, fobs, color='red', label='Picks', s=10)
    plt.show()
    plt.close()

    m, covm0, covm, f0_array, F_m = full_inversion(fobs, tobs, peaks_assos, mprior, num_iterations=4, sigma=5)
    #except:
    #    print('Error in full inversion for station:', sta, 'flight:', flight_num, 'date:', date)
    #    continue
    v0 = m[0]
    l = m[1]
    tprime0 = m[2]
    c = m[3]
    covm = np.sqrt(np.diag(covm))

    closest_index = np.argmin(np.abs(tprime0 - times))
    arrive_time = spec[:,closest_index]
    for i in range(len(arrive_time)):
        if arrive_time[i] < 0:
            arrive_time[i] = 0

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/with_c_quasi/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    qnum = plot_spectrogram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, covm, flight_num, middle_index, tarrive-start_time, closest_time, BASE_DIR, plot_show=False)
    qnum = "__"
    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/with_c_quasi/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, BASE_DIR)
    
    if rerun_fig == False:
        output.write(str(date)+','+str(flight_num)+','+str(sta)+','+str(closest_time)+','+str(v0)+','+str(l)+','+str(tprime0)+','+str(c)+','+str(f0_array)+','+str(covm)+','+str(qnum)+','+str(Tc)+','+str(c)+','+str(F_m)+',\n') 

    if rerun_fig == False:
        output.close()
