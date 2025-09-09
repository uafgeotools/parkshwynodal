import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
from src.doppler_funcs import *
from scipy.signal import find_peaks, spectrogram
import json
import obspy
from datetime import datetime, timezone
from pyproj import Proj
from plot_func import *
from obspy.clients.nrl import NRL

nrl = NRL()
equip = 'DH8A'
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

temp_correction = True #Flag use the temperature correction or not
rerun_fig = False #Flag rerun the figures without saving the inversion results = True

if rerun_fig == False:
    if temp_correction == True:
        output = open('output/' + equip + 'data_atmosphere_full.csv', 'a')
    else:
        output = open('output/' + equip + 'data_full.csv', 'a')
        
# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')


for li in file_in.readlines():

    text = li.split(',')
    flight_num = text[1]
    date = text[0]
    sta = text[9]
    time = float(text[5])

    equip_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_flights.csv'
    equip_data = pd.read_csv(equip_file, sep=",")
    equip_list = equip_data['equip']
    flight_list = equip_data['flight_id']


    for i_e in range(len(equip_list)):

        if str(flight_num) == str(flight_list[i_e]) and str(equip_list[i_e]) == str(equip):
            equip_go = True
            break
        else:
            equip_go = False
    if equip_go == False:
        continue

    ht = datetime.fromtimestamp(time, tz=timezone.utc)
    h = ht.hour

    alt = float(text[4])*0.0003048 #convert between feet and km
    x =  float(text[2]) 
    y = float(text[3])  

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)

    if temp_correction == True:
        input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'
        try:
            file =  open(input_files, 'r') #as file:
        except:
            print('No file for: ', date, flight_num, sta)
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
                    if abs(float(item['values'][i]) - float(alt)) < hold:
                        hold = abs(float(item['values'][i]) - float(alt))
                        z_index = i
        folder_spec = equip + '_spec_c'
        folder_spectrum = equip + '_spectrum_c'
        for item in data_list:
            if item['parameter'] == 'T':
                Tc = - 273.15 + float(item['values'][z_index])
    else:
        Tc = -2
        folder_spec =  equip + '_spec_cfc'
        folder_spectrum = equip + '_spectrum_cfc'
    c = speed_of_sound(Tc)
    sound_speed = c

    print(f"Speed of sound: {c} m/s")

    spec_dir = '/scratch/irseppi/nodal_data/plane_info/' + folder_spec +'/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'
    
    if os.path.exists(spec_dir) and rerun_fig == False:
        continue

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    time = flight_data['snapshot_id']
    timestamps = flight_data['snapshot_id']
    speed = flight_data['speed']
    altitude = flight_data['altitude']

    closest_x, closest_y, dist_km, closest_time, tarrive, alt, sp, elevation, speed_mps, height_m, dist_m, tmid = closest_approach_UTM(seismo_latitudes, seismo_longitudes, flight_latitudes, flight_longitudes, timestamps, altitude, speed, stations, elevations, c, sta)
    if closest_x == None:
        continue
            
    if equip == 'C185':
        #To set the initial window of arrival correct picks your start end Must use the tarrive time to get the correct data
        ta_old = calc_time(tmid,dist_m,height_m,343)
        start_time = ta_old - 120
        ht = datetime.fromtimestamp(ta_old, tz=timezone.utc)
    else:
        #Must use the tarrive time to get the correct data
        start_time = tarrive - 120
        ht = datetime.fromtimestamp(tarrive, tz=timezone.utc)
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
    except:
        continue 

    tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - 120, tr[2].stats.starttime + (mins * 60) + secs + 120)
    data = tr[2][:]
    fs = int(tr[2].stats.sampling_rate)
    title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
    torg = tr[2].times()

    # Compute spectrogram
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
    
    spec, MDF = remove_median(Sxx)
    
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    tprime0 = tarrive-start_time
    v0 = speed_mps
    l = np.sqrt(dist_m**2 + (height_m)**2)

    tf = np.arange(0, 240, 1)

    coords = doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time,make_picks=True) 

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
    if isinstance(sta, int):
        peaks = []
        p, _ = find_peaks(middle_column, distance = 7)
        corridor_width = (fs/2) / len(p) 
                        
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
        
        f0 = m[0]
        v0 = m[1]
        l = m[2]
        tprime0 = m[3]

    mprior = []
    mprior.append(v0)
    mprior.append(l)
    mprior.append(tprime0)       

    peaks, freqpeak =  overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time, tprime0, make_picks=True)
    f0_array = []
    w = len(peaks)
    for o in range(w):
        tprime = freqpeak[o]
        ft0p = peaks[o]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        mprior.append(f0)
        f0_array.append(f0)
    mprior = np.array(mprior)
                        
    corridor_width = 10

    peaks_assos = []
    fobs = []
    tobs = []

    for pp in range(len(peaks)):
        tprime = freqpeak[pp]
        ft0p = peaks[pp]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        ft = calc_ft(times,  tprime0, f0, v0, l, c)
        
        maxfreq = []
        coord_inv = []
        ttt = []

        f01 = f0 + corridor_width
        f02 = f0  - corridor_width
        upper = calc_ft(times,  tprime0, f01, v0, l, c)
        lower = calc_ft(times,  tprime0, f02, v0, l, c)

        for t_f in range(len(times)):

            try:      
                tt = spec[int(np.round(lower[t_f],0)):int(np.round(upper[t_f],0)), t_f]

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
                mtest = [f0,v0, l, tprime0]
                mtest,_, F_m = invert_f(mtest, coord_inv_array, c, num_iterations=4)
                ft = calc_ft(ttt,  mtest[3], mtest[0], mtest[1], mtest[2], c)
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


    if len(fobs) == 0:
        print('No picks for: ', date, flight_num, sta)
        continue

    tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, w, peaks_assos, make_picks=True)

    m, covm, f0_array, F_m = full_inversion(fobs, tobs, freqpeak, peaks, peaks_assos, tprime, tprime0, ft0p, v0, l, f0_array, mprior, c, w, 5)
    v0 = m[0]
    l = m[1]
    tprime0 = m[2]
    covm = np.sqrt(np.diag(covm))

    closest_index = np.argmin(np.abs(tprime0 - times))
    arrive_time = spec[:,closest_index]
    for i in range(len(arrive_time)):
        if arrive_time[i] < 0:
            arrive_time[i] = 0

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    qnum = plot_spectrgram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, covm, flight_num, middle_index, tarrive-start_time, closest_time, BASE_DIR, plot_show=False)

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, BASE_DIR)
    
    if rerun_fig == False:
        output.write(str(date)+','+str(flight_num)+','+str(sta)+','+str(closest_time)+','+str(tprime0)+','+str(v0)+','+str(l)+','+str(f0_array)+','+str(covm)+','+str(qnum)+','+str(Tc)+','+str(c)+','+str(F_m)+',\n') 

if rerun_fig == False:
    output.close()