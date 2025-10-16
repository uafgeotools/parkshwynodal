import numpy as np
import os
import gc
import concurrent.futures
from scipy.signal import spectrogram
from src.doppler_funcs import make_base_dir, load_waveform
from src.main_inv_fig_functions import remove_median, plot_spectrogram, plot_spectrum
import psutil

num_workers = os.cpu_count()

def inversion_process(line):
    text = line.split(',')
    date = text[0]
    month = int(date[4:6])
    day = date[6:8]
    flight_num = text[1]
    closest_time = float(text[5])
    sta = text[9]
    equip = text[10]

    input_file ='output/inv_results_ngt/' + equip + '_full_inv_results.txt'
    if os.path.exists(input_file) == False:
        return
    input_file = open(input_file, 'r')
    for inv_line in input_file.readlines():
        inv_text = inv_line.split(',')
        inv_flight_num = inv_text[1]
        inv_sta = inv_text[2]
        inv_closest_time = float(inv_text[3])
        if flight_num == inv_flight_num and sta == inv_sta and closest_time == inv_closest_time:
            v0 = float(inv_text[4])
            l = float(inv_text[5])
            t0 = float(inv_text[6])
            start_time = float(inv_text[7]) - t0 
            c = float(inv_text[8])
            f0_array = str(inv_text[9])
            f0_array = np.char.replace(f0_array, '[', '')
            f0_array = np.char.replace(f0_array, ']', '')
            f0_array = str(f0_array)
            f0_array = np.array(f0_array.split(' '))
            covm0 = inv_text[10]
            F_m = inv_text[13]
        else:
            continue
    input_file.close()

    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt_parallel/' + folder_spectrum + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+ '/'+str(sta)+'_' + str(closest_time) + '.png'
    if os.path.exists(DIR):
        return

    data, fs, t_wf, title = load_waveform(sta, start_time)
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant')

    spec, MDF = remove_median(Sxx)
    middle_index =  len(times) // 2

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt_parallel/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    _ = plot_spectrogram(data, fs, t_wf, title, spec, times, frequencies, t0, v0, l, c, f0_array, F_m, MDF, covm0, flight_num, middle_index, closest_time, BASE_DIR, plot_show=False, gt = False)

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt_parallel/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, times, frequencies, t0, l, c, f0_array, fs, closest_time, sta, BASE_DIR)

    # Explicitly delete large variables and collect garbage to free memory
    # Delete all variables and objects that may impact short-term memory
    del data, fs, t_wf, title
    del frequencies, times, Sxx, spec, MDF
    del v0, l, t0, start_time, c, f0_array, covm0, F_m
    del date, month, day, flight_num, closest_time, sta, equip
    del folder_spec, folder_spectrum, DIR, BASE_DIR

    gc.collect()
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")

# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

lines = file_in.readlines()
with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
   executor.map(inversion_process, lines)

