import numpy as np
import os
from datetime import datetime, timezone
from obspy.clients.nrl import NRL
from scipy.signal import spectrogram
from prelude import calc_time, make_base_dir, invert_f, full_inversion, get_speed_of_sound, get_sta_elevation, load_waveform
from main_inv_fig_functions import doppler_picks, overtone_picks, time_picks, remove_median, plot_spectrogram, plot_spectrum, get_auto_picks_full
jet = ['B737', 'B738', 'B739', 'B733', 'B763', 'B772', 'B77W', 'B788', 'B789', 'B744', 'B748', 'B77L', 'CRJ2', 'B732', 'A332', 'A359', 'E75S']

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
    day = date[6:8]
    flight_num = text[1]
    x =  float(text[2])  # UTM x-coordinate, meters
    y = float(text[3])  # UTM y-coordinate, meters
    dist_m = float(text[4])   # Distance in meters
    closest_time = float(text[5])
    alt = float(text[6]) 
    speed_mps = float(text[7])  # Speed in meters per second
    sta = text[9]
    equip = text[10]
    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    if mk_picks == False:
        file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/inversepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight_num) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight_num) + '.csv'
        if not os.path.exists(file_name):
            continue
    DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    if os.path.exists(DIR):
        continue

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
        error_file = open('output/inv_results/error_log.txt', 'a')
        error_file.write(f"Error loading waveform for {sta} on {date} at flight {flight_num}: {str(e)}\n")
        error_file.close()
        continue
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    coords, start_time = doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, tarrive, make_picks=mk_picks) 
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
    # Convert the list of coordinates to a numpy array
    coords_array = np.array(coords)

    # Find the index of the point in coords_array with the closest frequency to f0
    f0 = ((np.max(coords_array[:,1])+np.min(coords_array[:,1]))/2) - 20
    tprime0 = tarrive-start_time
    v0 = speed_mps
    height_m = alt - elev
    l = np.sqrt(dist_m**2 + (height_m)**2)

    mprior = []
    mprior.append(v0)
    mprior.append(l)
    mprior.append(tprime0)
    mprior.append(c)

    tf = np.arange(0, 240, 1)

    m0 = [f0, v0, l, tprime0, c]
    print('m0:', m0)
    sigma_prior = [20, 1, 1, 200, 1]
    m,_,_, F_m = invert_f(m0,sigma_prior, coords_array, num_iterations=8)
    tprime0 = m[3]

    sigma_f0 = 100
    sigma_v0 = 100
    sigma_l = 1000
    sigma_tprime0 = 200
    sigma_c = 100

    m0 = [f0, v0, l, tprime0, c]
    sigma_prior = [sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c]
    m,_,_, F_m = invert_f(m0,[sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c], coords_array, num_iterations=8)
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]
    mprior[2] = tprime0
    mprior[3] = c
    peaks, freqpeak =  overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time, tprime0, tarrive, make_picks=mk_picks)

    corridor_width = 8

    tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks,freqpeak, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax, equip)
    
    if len(fobs) == 0:
        continue

    for o in range(len(f0_array)):
        mprior.append(float(f0_array[o]))

    tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, len(peaks), peaks_assos, make_picks=mk_picks)

    if equip in jet:
        sigma_prior = [10, 10, 200, 10, 80]
    else:
        sigma_prior = [5, 10, 200, 30, 80]
    print('mprior:', mprior)
    m, covm0, covm, f0_array, F_m = full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations=4, sigma=3, off_diagonal=False)

    v0 = m[0]
    l = m[1]
    tprime0 = m[2]
    c = m[3]

    covm = np.sqrt(np.diag(covm))
    covm0 = np.sqrt(np.diag(covm0))

    closest_index = np.argmin(np.abs(tprime0 - times))
    arrive_time = spec[:,closest_index]
    for i in range(len(arrive_time)):
        if arrive_time[i] < 0:
            arrive_time[i] = 0

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    qnum = plot_spectrogram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, covm0, flight_num, middle_index, tarrive-start_time, closest_time, BASE_DIR, plot_show=False)
    qnum = "__"
    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, BASE_DIR)
    
    if rerun_fig == False:
        output = open('output/inv_results/' + equip + '_full_inv_results.csv', 'a')
        output.write(str(date)+','+str(flight_num)+','+str(sta)+','+str(closest_time)+','+str(v0)+','+str(l)+','+str(tprime0)+','+ str(start_time + tprime0) + ','+str(c)+','+str(f0_array)+','+str(covm0)+','+str(qnum)+','+str(Tc)+','+str(c)+','+str(F_m)+',\n') 
        output.close()
    else:
        continue  # Skip saving results if rerun_fig is True