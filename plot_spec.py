import numpy as np
import pandas as pd
import json
import obspy
import datetime
from datetime import datetime, timezone
from pyproj import Proj
from prelude import *
from scipy.signal import spectrogram
from plot_func import *
import os

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

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

    spec_dir = '/scratch/irseppi/nodal_data/plane_info/' + equip + '_spec_c/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.png'
    if os.path.exists(spec_dir):
        continue

    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + date + '_flights.csv', sep=",")
    flight = flight_data['flight_id']
    flight = flight.values.tolist()
    tailnumber = flight_data['aircraft_id']

    for i in range(len(stations)):
        if stations[i] == sta:
            elev = elevations[i]
            break
    
    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)

    input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'
    
    try:
        file =  open(input_files, 'r') 
    except:
        print('No tempurature file for: ', date, flight_num, sta)
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

    for item in data_list:
        if item['parameter'] == 'T':
            Tc = - 273.15 + float(item['values'][z_index])

    c = speed_of_sound(Tc)
    sound_speed = c
    tarrive = calc_time(closest_time,dist_m,(alt-elev),c) 


    #Must use the tarrive time to get the correct data
    ht = datetime.fromtimestamp(tarrive, tz=timezone.utc)

    wind = 120
    start_time = tarrive - wind
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

    try:
        p = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta)+".mseed"
        tr = obspy.read(p)
        tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - wind, tr[2].stats.starttime + (mins * 60) + secs + wind)
        data = tr[2][:]
        fs = int(tr[2].stats.sampling_rate)
        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
        torg = tr[2].times()
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
    except:
        try:
            p = "/scratch/irseppi/500sps/2019_0" + str(month) + "_" + str(day) + "/ZE_" + str(sta) + "_DPZ.msd"
            tr = obspy.read(p)

            tr.trim(tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs - wind, tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs + wind)

            data = tr[0][:]
            fs = int(tr[0].stats.sampling_rate)
            title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
            torg = tr[0].times()
        except:
            print(p)
            continue

    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
    try:
        spec, MDF = remove_median(Sxx)
    except:
        spec = Sxx
        MDF = np.full((0,250),0)
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    tf = np.arange(0, 240, 1)

    closest_index = np.argmin(np.abs(tarrive - times))
    arrive_time = spec[:,closest_index]
    for i in range(len(arrive_time)):
        if arrive_time[i] < 0:
            arrive_time[i] = 0

    # Plot settings and calculations
    vmin = np.min(arrive_time) #tprime0
    vmax = np.max(arrive_time)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

    ax1.plot(torg, data, 'k', linewidth=0.5)
    ax1.set_title(title)

    ax1.margins(x=0)
    ax1.set_position([0.125, 0.6, 0.775, 0.3])  # Move ax1 plot upwards

    # Plot spectrogram
    cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
    ax2.set_xlabel('Time (s)')

    ax2.axvline(x=tarrive-start_time, c = '#e41a1c', ls = '--',linewidth=0.5,label= r'$t_{i}$ = ' + "%.2f" % (tarrive-start_time) +' s')

    ax2.legend(loc='upper right',fontsize = 'small')
    ax2.set_ylabel('Frequency (Hz)')

    ax2.margins(x=0)
    ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

    plt.colorbar(mappable=cax, cax=ax3)
    ax3.set_ylabel('Relative Amplitude (dB)')

    ax2.margins(x=0)
    ax2.set_xlim(0, 240)
    ax2.set_ylim(0, int(fs/2))

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

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    print(BASE_DIR)
    make_base_dir(BASE_DIR)
    fig.savefig(BASE_DIR+'/'+str(closest_time)+'_'+str(flight_num)+'.png')
    plt.close()


    vmax = np.max(arrive_time)
    fig = plt.figure(figsize=(10,6))
    plt.grid()

    plt.plot(frequencies, spec[:,closest_index], c='#e41a1c')
    plt.xlim(0, int(fs/2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0,vmax*1.1)
    plt.xlabel('Frequency (Hz)', fontsize=17)
    plt.ylabel('Relative Amplitude at t = {:.2f} s (dB)'.format(tarrive-start_time), fontsize=17)

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    fig.savefig(BASE_DIR + '/'+str(sta)+'_' + str(closest_time) + '.png')
    plt.close()