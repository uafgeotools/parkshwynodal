import numpy as np
import pandas as pd
import obspy
import os
from pathlib import Path
from matplotlib import pyplot as plt
from datetime import datetime, timezone
from obspy.clients.nrl import NRL
from scipy.signal import spectrogram
from prelude import calc_ft, calc_time, speed_of_sound, calc_f0, make_base_dir, invert_f, full_inversion, get_speed_of_sound, get_sta_elevation, load_waveform
from main_inv_fig_functions import doppler_picks, overtone_picks, time_picks, remove_median, plot_spectrogram, plot_spectrum, get_auto_picks_1o, get_auto_picks_full
import shutil

nrl = NRL()
window = 120  # seconds before the arrival time to load the waveform
rerun_fig = True #Flag rerun the figures without saving the inversion results = True
mk_picks = False

# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

for li in file_in.readlines():
    text = li.split(',')
    date = text[0]
    month = int(date[4:6])
    day = int(date[6:8])
    flight_num = text[1]
    x =  float(text[2])  # UTM x-coordinate, meters
    y = float(text[3])  # UTM y-coordinate, meters
    dist_m = float(text[4])   # Distance in meters
    closest_time = float(text[5])
    alt = float(text[6]) 
    speed_mps = float(text[7])  # Speed in meters per second
    sta = text[9]
    equip = text[10]

    if rerun_fig == False:
        output = open('output/inv_results/' + equip + 'data_atmosphere_full.csv', 'a')

    elev = get_sta_elevation(sta)
    c, Tc = get_speed_of_sound(alt, closest_time, x, y)

    tarrive = calc_time(closest_time,dist_m,alt,c) 
    ht = datetime.fromtimestamp(tarrive, tz=timezone.utc)

    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    try:
        data, fs, torg, title = load_waveform(sta, (tarrive-window))
        frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant')
        spec, MDF = remove_median(Sxx)
    except Exception as e:
        print(Exception, e)
        continue
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    coords, start_time = doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, tarrive, make_picks=False) 
    coords_array = np.array(coords)
    if start_time is None or len(coords_array) == 0:
        continue

    elif equip == 'C185':
        start_time = start_time - 120

    if (tarrive - window) != start_time:
        data, fs, torg, title = load_waveform(sta, start_time)
        frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant')
        spec, MDF = remove_median(Sxx)

    if len(times) == 0 or len(frequencies) == 0 or len(Sxx) == 0:
        continue

    tprime0 = tarrive-start_time
    v0 = speed_mps
    height_m = alt - elev
    l = np.sqrt(dist_m**2 + (height_m)**2)

    tf = np.arange(0, 240, 1)

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

    coord_inv_array = get_auto_picks_1o(times, frequencies, spec, ft, corridor_width, m0, sigma_prior)
   
    m,_,_,F_m = invert_f(m0, sigma_prior, coord_inv_array, num_iterations=5, sigma=5)
        
    f0 = m[0]
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]

    peaks, freqpeak =  overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time, tprime0, tarrive, make_picks=True)

    corridor_width = (fs/2) / len(peaks) 
    if equip[0] == 'B' and equip[0:1] != 'BE':
        corridor_width = 3

    tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks,freqpeak, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax, equip)
    
    if len(fobs) == 0:
        continue

    tobs_hold = tobs.copy()

    tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, len(peaks), peaks_assos, make_picks=True)

    #if len(tobs) == len(tobs_hold):
    #    continue

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
