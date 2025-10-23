import os
import gc
import sys
import psutil
import numpy as np
from scipy.signal import spectrogram
from pathlib import Path

# Ensure repository root is on sys.path so local package 'src' can be imported
repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
from src.doppler_funcs import make_base_dir, invert_f, full_inversion, get_sta_elevation, load_waveform
from src.main_inv_fig_functions import time_picks, remove_median, plot_spectrogram, plot_spectrum, get_auto_picks_full

jet = ['B737', 'B738', 'B739', 'B733', 'B763', 'B772', 'B77W', 'B788', 'B789', 'B744', 'B748', 'B77L', 'CRJ2', 'B732', 'A332', 'A359', 'E75S']
rerun_fig = True #Flag rerun the figures without saving the inversion results = True
mk_picks = False

repo_path = '/home/irseppi/REPOSITORIES/parkshwynodal/'
fig_path = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt/'

# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open(repo_path + 'input/node_crossings_db_UTM.txt','r')

for line in file_in:
    text = line.split(',')
    date = text[0]
    month = int(date[4:6])
    day = date[6:8]
    flight_num = text[1]
    closest_time = float(text[5])
    sta = text[9]
    equip = text[10]
    alt = float(text[6]) 
    speed_gt = float(text[7]) 
    dist_m = float(text[4])   # Distance in meters
    elev = get_sta_elevation(sta)
    height_m = alt - elev
    distance_gt = np.sqrt(dist_m**2 + (height_m)**2) 

    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    DIR = fig_path + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    if os.path.exists(DIR):
        continue
    file_name = repo_path + 'input/data_picks/' + equip + '_data_picks/inversepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight_num) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight_num) + '.csv'
    if not os.path.exists(file_name):
        continue
    else:
        coords = []
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                coords.append((float(pick_data[0]), float(pick_data[1])))
            if len(pick_data) == 4:
                start_time = float(pick_data[2])
            else:
                file.close() 
                continue

        file.close()  

    coords_array = np.array(coords)
    if len(coords_array) == 0:
        continue

    elif equip == 'C185':
        start_time = start_time - 120

    c = 320 # Default speed of sound, average of dataset, m/s
    fa = np.max(coords_array[:, 1]) 
    fr = np.min(coords_array[:, 1])
    #insert method to get initial model here
    fm = (fa+fr)/2 

    #find the closest coordinate to f0
    closest_index = np.argmin(np.abs(coords_array[:, 1] - fm))
    f0 = coords_array[closest_index, 1] 
    t0 = coords_array[closest_index, 0]  
    t_hold = np.inf
    for i,t in enumerate(coords_array[:, 0]):
        if t != t0:
            if (t - t0) < t_hold:
                t_hold = abs(t - t0)
                second_index = i

    v0 = c*abs(fa-fr) / (2 * f0)
    slope = (coords_array[closest_index,1] - coords_array[second_index,1]) / (coords_array[closest_index,0] - coords_array[second_index,0])
    l = -((f0*v0**2/c)*(1-(v0/c)**2)**(-3/2))/slope 
    m0 = [f0, v0, l, t0, c]
    print('Initial model:', m0)

    data, fs, t_wf, title = load_waveform(sta, start_time)
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant')
    if len(times) == 0 or len(frequencies) == 0 or len(Sxx) == 0:
        continue

    spec, MDF = remove_median(Sxx)
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 
 
    m0 = [f0, v0, l, t0, c]
    sigma_prior = [40, 1, 1, 200, 1]
    m,_,_, F_m = invert_f(m0,sigma_prior, coords_array, num_iterations=3)
    m0[0] = m[0]
    m0[3] = m[3]

    tf = np.arange(0, 240, 1)

    sigma_f0 = 150
    sigma_v0 = 100
    sigma_l = 10000
    sigma_t0 = 200
    sigma_c = 100

    m0 = [f0, v0, l, t0, c]
    sigma_prior = [sigma_f0, sigma_v0, sigma_l, sigma_t0, sigma_c]
    m,_,_, F_m = invert_f(m0,[sigma_f0, sigma_v0, sigma_l, sigma_t0, sigma_c], coords_array, num_iterations=3)
    v0 = m[1]
    l = m[2]
    t0 = m[3]
    c = m[4]
    mprior = []
    mprior.append(v0)
    mprior.append(l)
    mprior.append(t0)
    mprior.append(c)


    output2 = repo_path + 'input/data_picks/' + equip + '_data_picks/overtonepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight_num) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight_num) + '.csv'
    if not os.path.exists(output2):
        continue
    else:
        peaks = []
        freqpeak = []
        with open(output2, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                peaks.append(float(pick_data[1]))
                freqpeak.append(float(pick_data[0]))
        file.close()  
    if len(peaks) <= 15:
        corridor_width = 10
    else:
        corridor_width = 5
    try:
        tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks,freqpeak, times, frequencies, spec, corridor_width, t0, v0, l, c, sigma_prior, vmax)
    except:
        continue

    if len(fobs) == 0:
        continue

    for o in range(len(f0_array)):
        mprior.append(float(f0_array[o]))

    tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, len(peaks), peaks_assos, make_picks=mk_picks)

    if abs(slope) < 1:
        sigma_prior = [10, 125, 15000, 30, 100]
    else:
        sigma_prior = [10, 30, 500, 30, 100]
    if equip in jet:
        sigma_prior = [100, 300, 50000, 100, 100]

    m, covm0, covm, f0_array, F_m = full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations=2, sigma=3, off_diagonal=False)

    v0 = m[0]
    l = m[1]
    t0 = m[2]
    c = m[3]

    covm = np.sqrt(np.diag(covm))
    covm0 = np.sqrt(np.diag(covm0))

    BASE_DIR = fig_path + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    qnum = plot_spectrogram(data, fs, t_wf, title, spec, times, frequencies, t0, v0, l, c, f0_array, F_m, MDF, covm0, flight_num, middle_index, closest_time, BASE_DIR, plot_show=False, gt = False)

    BASE_DIR = fig_path + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, times, frequencies, t0, l, c, f0_array, fs, closest_time, sta, BASE_DIR)

    if rerun_fig == False:
        output = open(repo_path + 'output/inv_results_ngt.txt', 'a')
        output.write(str(equip)+','+str(date)+','+str(flight_num)+','+str(sta)+','+str(closest_time)+','+str(v0)+','+str(l)+','+str(t0)+','+ str(start_time + t0) + ','+str(c)+','+str(f0_array)+','+str(covm0)+','+str(F_m)+',\n') 
        output.close()

    # Delete all variables and objects that may impact short-term memory
    del data, t_wf, frequencies, times, Sxx, spec, MDF
    del peaks, freqpeak, tobs, fobs, peaks_assos

    gc.collect()
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")