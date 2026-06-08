'''This script performs the full inversion for all non-jet aircraft in the 
denali nodal dataset, using the initial model calculated from the spectrograms. 
It also generates the spectrogram figures for each station. 
The results are saved to a text file and the figures are saved.
'''
import os
import sys
import psutil
import numpy as np
import concurrent.futures

from scipy.signal import spectrogram
from pathlib import Path


# --- Fix sys.path ---
repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.get_save_data import make_base_dir, load_waveform
from src.doppler_funcs import DopplerInversion as DI
from src.fig_func import SpecPlot, GetPicks


nw = os.cpu_count() # Number of workers for parallel processing
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

def get_start_time(month, day, flight_num, closest_time, sta, equip, repo_root):
    #Specific function to get picks acurately for Seppi 2025 paper
    picks_root = f'{repo_root}/input/data_picks/'
    base_dir = f'{picks_root}{equip}_data_picks/inversepicks/2019-0{month}-'
    base_dir += f'{day}/{flight_num}/{sta}/'
    file_name = f'{closest_time}_{flight_num}.csv'

    if Path(f'{base_dir}/{file_name}').exists():
        file = open(f'{base_dir}/{file_name}', 'r') 
        line = file.readline()
        pick_data = line.split(',')
        if len(pick_data) == 4:
            start_time = float(pick_data[2])
            if equip == 'C185':
                start_time = start_time - 120
        file.close()
        return start_time

def get_doppler_picks(
        spec_load, month, day, flight_num, closest_time, sta, equip, repo_root):
    
    picks_root = f'{repo_root}/input/data_picks/'
    base_dir = f'{picks_root}{equip}_data_picks/inversepicks/2019-0{month}-'
    base_dir += f'{day}/{flight_num}/{sta}/'
    file_name = f'{closest_time}_{flight_num}.csv'

    coords, _ = spec_load.single_doppler_data(file_name, base_dir)
    coords_array = np.array(coords)

    if len(coords_array) == 0:
        return
    
    return coords_array

def get_overtone_picks(
        spec_load, month, day, flight_num, closest_time, sta, equip, t0, v0, l, 
        c, sigma_prior, repo_root):

    picks_root = f'{repo_root}/input/data_picks/'
    base_dir = f'{picks_root}{equip}_data_picks/overtonepicks/2019-0{month}-'
    base_dir += f'{day}/{flight_num}/{sta}/'
    file_name = f'{closest_time}_{flight_num}.csv'

    peaks, freqpeak = spec_load.overtone_data(t0, file_name, base_dir)
   
    if len(peaks) <= 15:
        corridor_width = 10
    else:
        corridor_width = 5
    try:
        tobs, fobs, peaks_assos, fs_array = spec_load.auto_picks_full(
            peaks, freqpeak, corridor_width, t0, v0, l, c, sigma_prior)
    except:
        return

    if len(fobs) == 0:
        return
    
    DIR = f'{repo_root}/input/data_picks/{equip}_data_picks/timepicks/'
    DIR += f'2019-0{month}-{day}/{flight_num}/{sta}/'

    file_name = f'{closest_time}_{flight_num}.csv'
    tobs, fobs, peaks_assos = spec_load.final_data(tobs, fobs, fs_array, 
                                                   peaks_assos, file_name, DIR)

    return tobs, fobs, peaks_assos, fs_array

def inversion_process(
        month, day, flight_num, closest_time, sta, equip, repo_root):

    jet = ['B737', 'B738', 'B739', 'B733', 'B763', 'B772', 'B77W', 
       'B788', 'B789', 'B744', 'B748', 'B77L', 'CRJ2', 'B732', 
       'A332', 'A359', 'E75S'
       ]
    start_time = get_start_time(
        month, day, flight_num, closest_time, sta, equip, repo_root)

    data, sample_rate, t_wf, title = load_waveform(sta, start_time + 120, 120)
    
    # Compute spectrogram
    WIN_LEN = 1  # window length, in s
    NPER = int(WIN_LEN * sample_rate)
    frequencies, times, Sxx = spectrogram(
        data, sample_rate, scaling='density', nperseg=NPER,
        noverlap=int(NPER * .9), detrend='constant'
    )

    spec_plot = SpecPlot(Sxx, sample_rate, t_wf, data, times, frequencies, 
                    closest_time)
    
    spec, _ = spec_plot.remove_median()
    middle_index = len(times) // 2
    middle_column = spec[:, middle_index]
    vmin, vmax = 0, np.max(middle_column)
    
    spec_load = GetPicks(
        spec_plot, vmin, vmax, make_picks=False, save_picks=False)
    
    coords_array = get_doppler_picks(
        spec_load, month, day, flight_num, closest_time, sta, equip, repo_root)

    c = 320 # Default speed of sound, average of dataset, m/s
    fa = np.max(coords_array[:, 1])
    fr = np.min(coords_array[:, 1])
    #insert method to get initial model here
    fm = (fa+fr)/2 

    #find the closest coordinate to f0
    closest_index = np.argmin(np.abs(coords_array[:, 1] - fm))
    fs = coords_array[closest_index, 1] 
    t0 = coords_array[closest_index, 0]  
    t_hold = np.inf
    for i,t in enumerate(coords_array[:, 0]):
        if t != t0:
            if (t - t0) < t_hold:
                t_hold = abs(t - t0)
                second_index = i

    v = c*abs(fa-fr) / (2 * fs)
    y_diff = (coords_array[closest_index,1] - coords_array[second_index,1]) 
    xdiff = (coords_array[closest_index,0] - coords_array[second_index,0])

    slope = y_diff/xdiff
    d0 = -((fs*v**2/c)*(1-(v/c)**2)**(-3/2))/slope 
    m0 = [v, d0, t0, c, fs]
    print('Initial model:', m0)

    sigma_prior = [40, 1, 1, 200, 1]
    fobs = []
    tobs = []
    for t, f in coords_array:
        tobs.append(t)
        fobs.append(f)
    aircraft_inversion = DI(
        fobs, tobs, m0, sigma_prior, num_iterations=3, off_diagonal=False
        )
    # First inversion to refine model
    m, _, _, _, F_m = aircraft_inversion.full_inversion([len(fobs)])

    m0[4] = m[4]
    m0[2] = m[2]

    sigma_v = 100
    sigma_d0 = 10000
    sigma_t0 = 200
    sigma_c = 100
    sigma_fs = 150

    m0 = [v, d0, t0, c, fs]
    sigma_prior = [sigma_v, sigma_d0, sigma_t0, sigma_c, sigma_fs]
    aircraft_inversion = DI(
        fobs, tobs, m0, sigma_prior, num_iterations=3, off_diagonal=False)
    
    m, _, _, _, F_m = aircraft_inversion.full_inversion([len(fobs)])
    
    v = m[0]
    d0 = m[1]
    t0 = m[2]
    c = m[3]

    mprior = []
    mprior.append(v)
    mprior.append(d0)
    mprior.append(t0)
    mprior.append(c)

    tobs, fobs, peaks_assos, fs_array = get_overtone_picks(
        spec_load, month, day, flight_num, closest_time, sta, equip, t0, v, d0, 
        c, sigma_prior, repo_root)

    print('Overtone picks obtained for station', sta, 'and flight', flight_num)

    for o in range(len(fs_array)):
        mprior.append(float(fs_array[o]))

    if abs(slope) < 1:
        sigma_prior = [10, 125, 15000, 30, 100]
    else:
        sigma_prior = [10, 30, 500, 30, 100]
    if equip in jet:
        sigma_prior = [100, 300, 50000, 100, 100]

    print('Performing final inversion...')
    aircraft_inversion = DI(
        fobs, tobs, mprior, sigma_prior, num_iterations=4, off_diagonal=False)

    m, covm0, covm, fs_array, F_m = aircraft_inversion.full_inversion(
        peaks_assos, sigma=3)

    v = m[0]
    d0 = m[1]
    t0 = m[2]
    c = m[3]

    covm = np.sqrt(np.diag(covm))
    covm0 = np.sqrt(np.diag(covm0))

    return (
       spec_load, title, t0, v, d0, c, fs_array, F_m, covm0, flight_num, 
       middle_index
    )
    
def process_main(line, tracer, total_lines):
    print((tracer/total_lines)*100, '%')
    month, day, flight_num, closest_time, sta, equip = parse_line(line)
    fig_path = 'inversion_results_ngt/'
    folder_spec = equip + '_spec_c'

    DIR = f'{fig_path}{folder_spec}/2019-0{month}-{day}/{flight_num}/{sta}/'
    if os.path.exists(DIR):
        return

    (spec_load, title, t0, v, d0, c, fs_array, F_m, Cpost0, flight_num, 
     middle_index
     ) = inversion_process(
         month, day, flight_num, closest_time, sta, equip, repo_root)
    
    print('Inversion complete for station', sta, 'and flight', flight_num)

    BASE_DIR = f'{repo_root}/output/{fig_path}{folder_spec}'
    BASE_DIR += f'/20190{month}{day}/{flight_num}/{sta}/'

    make_base_dir(BASE_DIR)
    file_name = f'{BASE_DIR}{str(closest_time)}_{str(flight_num)}.png'
    
    spec_load.plot_spectrogram(title, t0, v, d0, c, fs_array, F_m, Cpost0, 
                               middle_index, file_name, plot_show=False, 
                               gt = False)

    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")

if __name__ == '__main__':
    # Loop through each station in text file that we already know 
    # comes within 2km of the nodes
    with open(f'{repo_root}/input/node_crossings_db_UTM.txt', 'r') as f:
        lines = f.readlines()

    tracer = [i for i in range(len(lines))]
    with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as executor:
        executor.map(process_main, lines, tracer, [len(lines)] * len(lines))
