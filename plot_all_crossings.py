import numpy as np
import os
import concurrent.futures
from datetime import datetime, timezone
from scipy.signal import spectrogram
import matplotlib.pyplot as plt
from src.doppler_funcs import make_base_dir, calc_time, get_speed_of_sound, get_sta_elevation, load_waveform
from src.main_inv_fig_functions import remove_median
from matplotlib.ticker import MaxNLocator

num_workers = os.cpu_count()


def plot_spectrogram(data, fs, torg, title, spec, times, frequencies, arrive_time, MDF, flight, middle_index, closest_time, dir_name):

    vmin = 0 
    vmax = np.max(spec)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     
    ax1.plot(torg, data, 'k', linewidth=0.5)
    ax1.set_title(title)

    ax1.margins(x=0)
    ax1.set_position([0.125, 0.6, 0.775, 0.3]) 
    # Plot spectrogram
    cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
    ax2.set_xlabel('Time (s)')

    ax2.set_ylabel('Frequency (Hz)')

    ax2.margins(x=0)
    ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

    # Set colorbar with integer ticks only
    cbar = plt.colorbar(mappable=cax, cax=ax3)
    cbar.locator = MaxNLocator(integer=True)
    cbar.update_ticks()
    ax3.set_ylabel('Relative Amplitude (dB)')

    ax2.margins(x=0)
    ax2.set_xlim(0, 240)
    ax2.set_ylim(0, int(fs/2))

    ax1.tick_params(axis='both', which='major', labelsize=9)
    ax2.tick_params(axis='both', which='major', labelsize=9)
    ax3.tick_params(axis='both', which='major', labelsize=9)

    # Plot overlay
    spec2 = 10 * np.log10(MDF)
    middle_column2 = spec2[:, middle_index]
    vmin2 = np.min(middle_column2)
    vmax2 = np.max(middle_column2)

    # Create ax4 and plot on the same y-axis as ax2
    ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) 
    ax4.plot(middle_column2, frequencies, c='#ff7f00')  
    ax4.set_ylim(0, int(fs/2))
    ax4.set_xlim(vmax2*1.1, vmin2) 
    ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    ax4.grid(axis='y')

    fig.savefig(dir_name+'/'+str(closest_time)+'_'+str(flight)+'.png')
    plt.close()

def plot_spectrum(spec, frequencies, fs, closest_index, closest_time, sta, dir_name):

    vmax = np.max(spec[:,closest_index])
    fig = plt.figure(figsize=(10,6))
    plt.grid()

    plt.plot(frequencies, spec[:,closest_index], c='#377eb8')
  
    plt.xlim(0, int(fs/2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0,vmax*1.1)
    plt.xlabel('Frequency (Hz)', fontsize=17)
    plt.ylabel('Relative Amplitude at t = 120 s (dB)', fontsize=17)

    fig.savefig(dir_name + '/'+str(sta)+'_' + str(closest_time) + '.pdf')
    plt.close()
def load_plot_spectrogram(sta, date, flight_num, tarrive, closest_time):
    window = 120  # seconds before the arrival time to load the waveform
    try:
        data, fs, torg, title = load_waveform(sta, (tarrive-window))
        frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend='constant')
        spec, MDF = remove_median(Sxx)
    except Exception as e:
        error_file = open('output/spec_error)_log.txt', 'a')
        error_file.write(f"Error loading waveform for {sta} on {date} at flight {flight_num}: {str(e)}\n")
        error_file.close()

    middle_index =  len(times) // 2
    base_dir = '/scratch/irseppi/nodal_data/plane_info/spec_no_inv/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)
    make_base_dir(base_dir)
    plot_spectrogram(data, fs, torg, title, spec, times, frequencies, 120, MDF, flight_num, middle_index, closest_time, base_dir)
    
    BASE_DIR =  '/scratch/irseppi/nodal_data/plane_info/spec_no_inv/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, frequencies, fs, middle_index, closest_time, sta, BASE_DIR)
    return BASE_DIR
# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
sta_list = []
date_list = []
flight_list = []
tarrive_list = []
closest_time_list = []
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

    DIR = '/scratch/irseppi/nodal_data/plane_info/spec_no_inv/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.png'
    if os.path.exists(DIR):
        if os.path.isfile(DIR):
            continue

    elev = get_sta_elevation(sta)
    c, Tc = get_speed_of_sound(alt, closest_time, x, y)

    tarrive = calc_time(closest_time,dist_m,alt,c) 
    sta_list.append(sta)
    date_list.append(date)
    flight_list.append(flight_num)
    tarrive_list.append(tarrive)
    closest_time_list.append(closest_time)

with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
    dir = executor.map(load_plot_spectrogram, sta_list, date_list, flight_list, tarrive_list, closest_time_list)