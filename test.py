import numpy as np
import os
import gc
from obspy.clients.nrl import NRL
from scipy.signal import spectrogram
from src.doppler_funcs import make_base_dir, invert_f, full_inversion, get_sta_elevation, load_waveform
from src.main_inv_fig_functions import time_picks, remove_median, plot_spectrogram, plot_spectrum, get_auto_picks_full
import psutil

paper_figures = ['C185_20190221_529754214_1550781331.5739982_1011_C185', 'B190_20190227_530696852_1551228121.0402486_1049_B190', 'B737_20190225_530339730_1551061570.9016998_1126_B737', 'B737_20190304_531697514_1551714047.0320563_1122_B737', 'B737_20190304_531711629_1551719807.3910785_1072_B737','B763_20190214_528407493_1550165581.4383187_1284_B763','C46_20190222_529805251_1550803683.768247_1007_C46', 'C185_20190221_529754214_1550777713.1677284_1020_C185', 'DH8A_20190214_528445164_1550158750.7401662_1173_DH8A', 'R44_20190213_528293430_1550089022.9259548_1007_R44']
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
    closest_time = float(text[5])
    sta = text[9]
    equip = text[10]

    input_file ='output/inv_results_ngt/' + equip + '_full_inv_results.txt'
    if os.path.exists(input_file) == False:
        continue
    else:
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

    file_check = str(equip)+'_'+ str(date) +'_'+str(flight_num)+'_' + str(closest_time) + '_' + str(sta) + '_' + str(equip)
    if file_check not in paper_figures:
        continue
    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt_tttest/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    #if os.path.exists(DIR):
    #    continue

    if equip == 'C185':
        start_time = start_time - 120

    data, fs, t_wf, title = load_waveform(sta, start_time)
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant')
  
    if len(times) == 0 or len(frequencies) == 0 or len(Sxx) == 0:
        continue
    spec, MDF = remove_median(Sxx)
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt_tttest/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    qnum = plot_spectrogram(data, fs, t_wf, title, spec, times, frequencies, t0, v0, l, c, f0_array, F_m, MDF, covm0, flight_num, middle_index, closest_time, BASE_DIR, plot_show=False, gt = False)
    qnum = "__"
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage spec 1: {mem:.2f} MB")

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt_tttest/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, times, frequencies, t0, l, c, f0_array, fs, closest_time, sta, BASE_DIR)
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage spec 2: {mem:.2f} MB")

    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage post: {mem:.2f} MB")

    del data, fs, t_wf, title
    del frequencies, times, Sxx, spec, MDF
    del covm0, f0_array, F_m, BASE_DIR
    del date, month, day, flight_num, closest_time, sta, equip
    del folder_spec, folder_spectrum, DIR
    del start_time, c, t0, v0, l

    gc.collect()
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")