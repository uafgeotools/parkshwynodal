'''This script performs the full inversion for all non-jet aircraft in the 
denali nodal dataset, using the initial model calculated from the spectrograms. 
It also generates the spectrogram figures for each station. 
The results are saved to a text file and the figures are saved.
'''
import os
import sys
import psutil
import numpy as np

from scipy.signal import spectrogram
from pathlib import Path

import concurrent.futures

# --- Fix sys.path ---
repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.doppler_funcs import (
    make_base_dir, invert_f, full_inversion, load_waveform)

from src.main_inv_fig_functions import (
    time_picks, remove_median, plot_spectrogram, 
    get_auto_picks_full)

num_workers = os.cpu_count()
picks_root = f'{repo_root}/input/data_picks/'

# Loop through each station in text file that we already 
# know comes within 2km of the nodes
file_in = open(f'{repo_root}/input/node_crossings_db_UTM.txt','r')

def parse_line(line):
    text = line.split(',')
    date = text[0]
    month = int(date[4:6])
    day = date[6:8]
    flight_num = text[1]
    closest_time = float(text[5])
    sta = text[9]
    equip = text[10]
    return month, day, flight_num, closest_time, sta, equip

def get_doppler_picks(month, day, flight_num, closest_time, sta, equip, repo_root):
    picks_root = f'{repo_root}/input/data_picks/'
    file_name = f'{picks_root}{equip}_data_picks/inversepicks/2019-0{month}-{day}'
    file_name += f'/{flight_num}/{sta}/{closest_time}_{flight_num}.csv'

    if not os.path.exists(file_name):
        return

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
                return

        file.close()
    
    coords_array = np.array(coords)
    if len(coords_array) == 0:
        return
    
    elif equip == 'C185':
        start_time = start_time - 120

    return coords_array, start_time 

def get_overtone_picks(
        month, day, flight_num, closest_time, sta, equip, start_time, times, 
        frequencies, spec, t0, v0, l, c, sigma_prior,vmax, repo_root):

    picks_root = f'{repo_root}/input/data_picks/'
    output2 = f'{picks_root}{equip}_data_picks/overtonepicks/2019-0{month}-{day}'
    output2 += f'/{flight_num}/{sta}/{closest_time}_{flight_num}.csv'

    if not os.path.exists(output2):
        return
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
        tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(
            peaks, freqpeak, times, frequencies, spec, corridor_width, 
            t0, v0, l, c, sigma_prior,vmax)
    except:
        return

    if len(fobs) == 0:
        return

    tobs, fobs, peaks_assos = time_picks(
        month, day, flight_num, sta, equip, tobs, fobs, closest_time, 
        start_time, spec, times, frequencies, 0, vmax, len(peaks), 
        peaks_assos, make_picks=False)

    return tobs, fobs, peaks_assos, f0_array

def inversion_process(month, day, flight_num, closest_time, sta, equip, repo_root):

    jet = ['B737', 'B738', 'B739', 'B733', 'B763', 'B772', 'B77W', 
       'B788', 'B789', 'B744', 'B748', 'B77L', 'CRJ2', 'B732', 
       'A332', 'A359', 'E75S'
       ]
    
    coords_array, start_time = get_doppler_picks(
        month, day, flight_num, closest_time, sta, equip, repo_root)

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
    y_diff = (coords_array[closest_index,1] - coords_array[second_index,1]) 
    xdiff = (coords_array[closest_index,0] - coords_array[second_index,0])

    slope = y_diff/xdiff
    l = -((f0*v0**2/c)*(1-(v0/c)**2)**(-3/2))/slope 
    m0 = [f0, v0, l, t0, c]
    print('Initial model:', m0)

    data, fs, t_wf, title = load_waveform(sta, (start_time+120))
    frequencies, times, Sxx = spectrogram(
        data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, 
        detrend = 'constant')
    
    if len(times) == 0 or len(frequencies) == 0 or len(Sxx) == 0:
        return

    spec, MDF = remove_median(Sxx)
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]

    vmax = np.max(middle_column) 
 
    m0 = [f0, v0, l, t0, c]
    sigma_prior = [40, 1, 1, 200, 1]
    m,_,_, F_m = invert_f(m0, sigma_prior, coords_array, num_iterations=3)
    m0[0] = m[0]
    m0[3] = m[3]

    sigma_f0 = 150
    sigma_v0 = 100
    sigma_l = 10000
    sigma_t0 = 200
    sigma_c = 100

    m0 = [f0, v0, l, t0, c]
    sigma_prior = [sigma_f0, sigma_v0, sigma_l, sigma_t0, sigma_c]
    m,_,_, F_m = invert_f(
        m0,[sigma_f0, sigma_v0, sigma_l, sigma_t0, sigma_c],
        coords_array, num_iterations=3)
    v0 = m[1]
    l = m[2]
    t0 = m[3]
    c = m[4]
    mprior = []
    mprior.append(v0)
    mprior.append(l)
    mprior.append(t0)
    mprior.append(c)

    tobs, fobs, peaks_assos, f0_array = get_overtone_picks(
        month, day, flight_num, closest_time, sta, equip, start_time, times, 
        frequencies, spec, t0, v0, l, c, sigma_prior,vmax, repo_root)


    for o in range(len(f0_array)):
        mprior.append(float(f0_array[o]))

    if abs(slope) < 1:
        sigma_prior = [10, 125, 15000, 30, 100]
    else:
        sigma_prior = [10, 30, 500, 30, 100]
    if equip in jet:
        sigma_prior = [100, 300, 50000, 100, 100]

    
    m, covm0, covm, f0_array, F_m = full_inversion(
        fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations=2, sigma=3, 
        off_diagonal=False)


    v0 = m[0]
    l = m[1]
    t0 = m[2]
    c = m[3]

    covm = np.sqrt(np.diag(covm))
    covm0 = np.sqrt(np.diag(covm0))

    return (
       data, fs, t_wf, title, spec, times, frequencies, t0, v0, l, c, f0_array, 
        F_m, MDF, covm0, flight_num, middle_index, closest_time
    )
    
def process_main(line, tracer, total_lines):
    print((tracer/total_lines)*100, '%')
    month, day, flight_num, closest_time, sta, equip = parse_line(line)
    fig_path = 'inversion_results_ngt/'
    folder_spec = equip + '_spec_c'

    DIR = f'{fig_path}{folder_spec}/2019-0{month}-{day}/{flight_num}/{sta}/'
    if os.path.exists(DIR):
        return
    (data, fs, t_wf, title, spec, times, frequencies, t0, v0, l, c, f0_array, 
     F_m, MDF, covm0, flight_num, middle_index, closest_time
     ) = inversion_process(month, day, flight_num, closest_time, sta, equip, repo_root)

    BASE_DIR = f'{repo_root}/output/{fig_path}{folder_spec}'
    BASE_DIR += f'/20190{month}{day}/{flight_num}/{sta}/'

    make_base_dir(BASE_DIR)
    plot_spectrogram(
        data, fs, t_wf, title, spec, times, frequencies, t0, v0, l, c, f0_array, 
        F_m, MDF, covm0, flight_num, middle_index, closest_time, BASE_DIR, 
        plot_show=False, gt = False)

    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")
if __name__ == '__main__':
    # Loop through each station in text file that we already know 
    # comes within 2km of the nodes
    with open(f'{repo_root}/input/node_crossings_db_UTM.txt', 'r') as f:
        lines = f.readlines()

    tracer = [i for i in range(len(lines))]
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        executor.map(process_main, lines, tracer, [len(lines)] * len(lines))

