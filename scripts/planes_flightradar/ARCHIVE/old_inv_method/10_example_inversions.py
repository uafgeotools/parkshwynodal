import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime

from src.doppler_funcs import *
from scipy.signal import find_peaks, spectrogram
from pathlib import Path
show_process = False
auto_peak_pick = False
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
elevations = seismo_data['Elevation']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430,527937367,529741194,529776675,529179112,530165646,531605202,531715679,529805251,529948401] 
time = [1551066051,1550172833,1550168070,1550165577,1550089044,1549912188,1550773710,1550787637,1550511447,1550974151,1551662362,1551736354,1550803701,1550867033] 
sta = [1022,1272,1173,1283,1004,"CCB","F6TP","F4TN","F3TN","F7TV",1010,1021,1006,1109] 
day = [25,14,14,14,13,11,21,21,18,24,4,4,22,22]
month = [2,2,2,2,2,2,2,2,2,2,3,3,2,2]

for n in range(0,13):
    ht = datetime.utcfromtimestamp(time[n])
    mins = ht.minute
    secs = ht.second
    h = ht.hour
    tim = 120	
    h_u = str(h+1)
    if h < 23:			
        day2 = str(day[n])
        if h < 10:
            h_u = '0'+str(h+1)
            h = '0'+str(h)
        else:
            h_u = str(h+1)
            h = str(h)
    else:
        h_u = '00'
        day2 = str(day[n]+1)
    if len(str(day[n])) == 1:
        day[n] = '0'+str(day[n])
        day2 = day[n]
    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/20190'+str(month[n])+str(day[n])+'_positions/20190'+str(month[n])+str(day[n])+'_'+str(flight_num[n])+'.csv', sep=",")

    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    tm = flight_data['snapshot_id']
    speed = flight_data['speed']
    alt = flight_data['altitude']
    head = flight_data['heading']
    for line in range(len(tm)):
        if str(tm[line]) == str(time[n]):
            speed = flight_data['speed'][line]
            speed_mps = speed * 0.514444
            alt = flight_data['altitude'][line]
            alt_m = alt * 0.3048
            for y in range(len(station)):
                if str(station[y]) == str(sta[n]):
                    dist = distance(seismo_latitudes[y], seismo_longitudes[y], flight_latitudes[line], flight_longitudes[line])	
                    if isinstance(sta[n], str):
                        day_of_year = str((ht - datetime(2019, 1, 1)).days + 1)
	
                        p = "/aec/wf/2019/0"+day_of_year+"/"+str(sta[n])+".*Z.20190"+day_of_year+"000000+"
                        tr = obspy.read(p)
                        tr[0].trim(tr[0].stats.starttime +(int(h) *60 *60) + (mins * 60) + secs - tim, tr[0].stats.starttime +(int(h) *60 *60) + (mins * 60) + secs + tim)
                        data = tr[0][:]
                        fs = int(tr[0].stats.sampling_rate)
                        title    = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'						
                        t_wf                  = tr[0].times()
                    else:
                        p = "/scratch/naalexeev/NODAL/2019-0"+str(month[n])+"-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-0"+str(month[n])+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(station[y])+".mseed"
                        tr = obspy.read(p)
                        tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
                        data = tr[2][:]
                        fs = int(tr[2].stats.sampling_rate)
                        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
                        t_wf = tr[2].times()
                      
                    clat, clon, dist_m, tmid = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
                    dist_km = dist_m / 1000
                    tarrive = tim + (time[n] - calc_time(tmid,dist_m,alt_m))
                    
                    # Compute spectrogram
                    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
                    
                    a, b = Sxx.shape

                    MDF = np.zeros((a,b))
                    for row in range(len(Sxx)):
                        m = len(Sxx[row])
                        p = sorted(Sxx[row])
                        median = p[int(m/2)]
                        for col in range(m):
                            MDF[row][col] = median
                    spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

                    ty = False
                    if ty == True:
                        if isinstance(sta[n], int):
                            spec = np.zeros((a,b))
                            for col in range(0,b):
                                p = sorted(Sxx[:, col])
                                median = p[int(len(p)/2)]

                                for row in range(len(Sxx)):
                                    spec[row][col] = 10 * np.log10(Sxx[row][col]) - ((10 * np.log10(MDF[row][col])) + ((10*np.log10(median))))

                    middle_index = len(times) // 2
                    middle_column = spec[:, middle_index]
                    vmin = 0  
                    vmax = np.max(middle_column) 
                    p, _ = find_peaks(middle_column, distance=10)

                    output1 = '/scratch/irseppi/nodal_data/plane_info/inversepicks/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.csv'
                    
                    if Path(output1).exists():
                        plt.figure()
                        plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                        pick_data = pd.read_csv(output1, header=None)

                        plt.scatter(pick_data.iloc[:, 0], pick_data.iloc[:, 1], color='black', marker='x')
                        plt.show()


                        con = input("Do you want to overwrite your previously stored data with your picks? (y or n)")
                        while con != 'y' and con != 'n':
                            print('Invalid input. Please enter y or n.')
                            con = input("Do you want to overwrite your previously stored data with your picks? (y or n)")
                        if con == 'y':
                            r1 = open(output1,'w')
                            coords = []
                            plt.figure()
                            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

                            def onclick(event):
                                global coords
                                coords.append((event.xdata, event.ydata))
                                plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                                plt.draw() 
                                print('Clicked:', event.xdata, event.ydata)  
                                r1.write(str(event.xdata) + ',' + str(event.ydata) + ',\n')

                            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                            plt.show(block=True)
                            r1.close()

                            pick_again = input("Do you want to repick you points? (y or n)")
                            while pick_again == 'y':
                                r1 = open(output1,'w')
                                coords = []
                                plt.figure()
                                plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)  
                                def onclick(event):
                                    global coords
                                    coords.append((event.xdata, event.ydata))
                                    plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                                    plt.draw() 
                                    print('Clicked:', event.xdata, event.ydata)  
                                    r1.write(str(event.xdata) + ',' + str(event.ydata) + ',\n') 
                                cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                                plt.show(block=True)
                                r1.close()
                                pick_again = input("Do you want to repick you points? (y or n)")
                        else:
                            coords = []
                            with open(output1, 'r') as file:
                                for line in file:
                                    # Split the line using commas
                                    pick_data = line.split(',')
                                    coords.append((float(pick_data[0]), float(pick_data[1])))
                            file.close()  # Close the file after reading

                    else:
                        BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversepicks/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'
                        make_base_dir(BASE_DIR)
                        r1 = open(output1,'w')
                        coords = []
                        plt.figure()
                        plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

                        def onclick(event):
                            global coords
                            coords.append((event.xdata, event.ydata))
                            plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                            plt.draw() 
                            print('Clicked:', event.xdata, event.ydata)  
                            r1.write(str(event.xdata) + ',' + str(event.ydata) + ',\n')

                        cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                        plt.show(block=True)
                        r1.close()

                        pick_again = input("Do you want to repick you points? (y or n)")
                        while pick_again == 'y':
                            r1 = open(output1,'w')
                            coords = []
                            plt.figure()
                            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                            def onclick(event):
                                global coords
                                coords.append((event.xdata, event.ydata))
                                plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                                plt.draw() 
                                print('Clicked:', event.xdata, event.ydata)  
                                r1.write(str(event.xdata) + ',' + str(event.ydata) + ',\n')
                            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                            plt.show(block=True)
                            r1.close()
                            pick_again = input("Do you want to repick you points? (y or n)")
                    # Convert the list of coordinates to a numpy array
                    coords_array = np.array(coords)

                    if n == 0:
                        t0 = 112
                        f0 = 115
                        v0 = 68
                        l = 2135

                    if n == 1:
                        f0 = 110
                        t0 = 107
                        v0 = 100
                        l = 2700

                    if n == 2:
                        f0 = 131
                        t0 = 93
                        v0 = 139
                        l = 4650

                    if n == 3:
                        f0 = 121
                        t0 = 116
                        v0 = 142
                        l = 2450

                    if n == 4:
                        f0 = 120
                        t0 = 140
                        v0 = 64
                        l = 580

                    if n == 5:
                        f0 = 17.5
                        t0 = 123
                        v0 = 112
                        l = 1150

                    if n == 6:
                        f0 = 36
                        t0 = 133
                        v0 = 92
                        l = 2400

                    if n == 7:
                        f0 = 26		
                        t0 = 122
                        v0 = 126
                        l = 3000

                    if n == 8:
                        f0 = 87.7
                        t0 = 100
                        v0 = 67
                        l = 2300

                    if n == 9:
                        f0 = 26
                        t0 = 114
                        v0 = 144
                        l = 1900
                    if n > 9:

                        f0 = fs/4
                        t0 = tarrive
                        v0 = speed_mps
                        l = np.sqrt(dist_m**2 + (alt_m-elevations[y])**2)

                    c = 343
                    m0 = [f0, v0, l, t0]

                    m,covm = invert_f(m0, coords_array, num_iterations=8)
                    f0 = m[0]
                    v0 = m[1]
                    l = m[2]
                    t0 = m[3]
                    
                    ft = calc_ft(times, t0, f0, v0, l, c)
                    if isinstance(sta[n], int):
                        peaks = []
                        p, _ = find_peaks(middle_column, distance = 7)
                        corridor_width = (fs/2) / len(p) 
                                        
                        if len(p) == 0:
                            corridor_width = fs/4

                        coord_inv = []
                        if show_process == True:
                            plt.figure()
                            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

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
                            if show_process == True:
                                plt.scatter(times[t_f], max_amplitude_frequency, color='black', marker='x')
                                plt.scatter(times[t_f], upper, color='pink', marker='x')
                                plt.scatter(times[t_f], lower, color='pink', marker='x')
                        if show_process == True:
                            plt.show()

                        coord_inv_array = np.array(coord_inv)

                        m,_ = invert_f(m0, coord_inv_array, num_iterations=12)
                        f0 = m[0]
                        v0 = m[1]
                        l = m[2]
                        t0 = m[3]

                        ft = calc_ft(times, t0, f0, v0, l, c)
                        
                        delf = np.array(ft) - np.array(peaks)
                        
                        new_coord_inv_array = []
                        for i in range(len(delf)):
                            if np.abs(delf[i]) <= 3:
                                new_coord_inv_array.append(coord_inv_array[i])
                        coord_inv_array = np.array(new_coord_inv_array)

                        if show_process == True:
                            plt.figure()
                            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                            plt.scatter(coord_inv_array[:,0], coord_inv_array[:,1], color='black', marker='x')
                            plt.show()

                        m,covm = invert_f(m0, coord_inv_array, num_iterations=12, sigma=5)
                        
                        f0 = m[0]
                        v0 = m[1]
                        l = m[2]
                        t0 = m[3]

                        ft = calc_ft(times, t0, f0, v0, l, c)
                        if show_process == True:
                            plt.figure()
                            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                            plt.plot(times, ft, color='g')
                            plt.show()        

                    if auto_peak_pick == True:
                        # Find the closest time value to m[3]
                        closest_time_index = np.argmin(np.abs(times - m[3]))

                        # Extract the corresponding column from the spectrogram
                        col = spec[:, closest_time_index]
                        peaks, _ = find_peaks(col,prominence=15) #distance = 10) 

                        plt.figure()
                        plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                        plt.axvline(x=t0, c = 'g', ls = '--')
                        plt.plot(col)
                        for p in peaks:
                            plt.scatter(times[closest_time_index], frequencies[p], color='black', marker='x')
                            plt.scatter(p, col[p], color='black', marker='x')
                        plt.show()
                    else:
                        output2 = '/scratch/irseppi/nodal_data/plane_info/overtonepicks/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.csv'
                        if Path(output2).exists():
                            plt.figure()
                            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                            plt.axvline(x=t0, c = '#377eb8', ls = '--')
                            pick_data = pd.read_csv(output2, header=None)

                            plt.scatter(pick_data.iloc[:, 0], pick_data.iloc[:, 1], color='black', marker='x')
                            plt.show()

                            con = input("Do you want to overwrite your previously stored overtone picks? (y or n)")
                            while con != 'y' and con != 'n':
                                print('Invalid input. Please enter y or n.')
                                con = input("Do you want to overwrite your previously stored overtone picks? (y or n)")
                            if con == 'y':
                                r2 = open(output2,'w')
                                peaks = []
                                freqpeak = []
                                plt.figure()
                                plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                                plt.axvline(x=t0, c = '#377eb8', ls = '--')
                                def onclick(event):
                                    global coords
                                    peaks.append(event.ydata)
                                    freqpeak.append(event.xdata)
                                    plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                                    plt.draw() 
                                    print('Clicked:', event.xdata, event.ydata)  
                                    r2.write(str(event.xdata) + ',' + str(event.ydata) + ',\n')

                                cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                                plt.show(block=True)
                                r2.close()

                                pick_again = input("Do you want to repick you points? (y or n)")
                                while pick_again == 'y':
                                    r2 = open(output2,'w')
                                    peaks = []
                                    freqpeak = []
                                    plt.figure()
                                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                                    plt.axvline(x=t0, c = '#377eb8', ls = '--')
                                    def onclick(event):
                                        global coords
                                        peaks.append(event.ydata)
                                        freqpeak.append(event.xdata)
                                        plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                                        plt.draw() 
                                        print('Clicked:', event.xdata, event.ydata)  
                                        r2.write(str(event.xdata) + ',' + str(event.ydata) + ',\n')
                                    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                                    plt.show(block=True)
                                    r2.close()
                                    pick_again = input("Do you want to repick you points? (y or n)")

                            elif con == 'n':
                                peaks = []
                                freqpeak = []
                                with open(output2, 'r') as file:
                                    for line in file:
                                        # Split the line using commas
                                        pick_data = line.split(',')
                                        peaks.append(float(pick_data[1]))
                                        freqpeak.append(float(pick_data[0]))
                                file.close()  # Close the file after reading

                        else:
                            BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/overtonepicks/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'
                            make_base_dir(BASE_DIR)
                            r2 = open(output2,'w')
                            peaks = []
                            freqpeak = []
                            plt.figure()
                            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                            plt.axvline(x=t0, c = '#377eb8', ls = '--')

                            def onclick(event):
                                global coords
                                peaks.append(event.ydata)
                                freqpeak.append(event.xdata)
                                plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                                plt.draw() 
                                print('Clicked:', event.xdata, event.ydata)  
                                r2.write(str(event.xdata) + ',' + str(event.ydata) + ',\n')
                            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
                            plt.show(block=True)
                            r2.close()
                            pick_again = input("Do you want to repick you points? (y or n)")
                            while pick_again == 'y':
                                r2 = open(output2,'w')
                                peaks = []
                                freqpeak = []
                                plt.figure()
                                plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                                plt.axvline(x=t0, c = '#377eb8', ls = '--')
                                def onclick(event):
                                    global coords
                                    peaks.append(event.ydata)
                                    freqpeak.append(event.xdata)
                                    plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                                    plt.draw() 
                                    print('Clicked:', event.xdata, event.ydata)  
                                    r2.write(str(event.xdata) + ',' + str(event.ydata) + ',\n')
                                cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                                plt.show(block=True)
                                r2.close()
                                pick_again = input("Do you want to repick you points? (y or n)")
                    closest_index = np.argmin(np.abs(t0 - times))
                    arrive_time = spec[:,closest_index]
                    for i in range(len(arrive_time)):
                        if arrive_time[i] < 0:
                            arrive_time[i] = 0
                    vmin = np.min(arrive_time) 
                    vmax = np.max(arrive_time) 

                    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

                    ax1.plot(t_wf, data, 'k', linewidth=0.5)
                    ax1.set_title(title)

                    ax1.margins(x=0)
                    ax1.set_position([0.125, 0.6, 0.775, 0.3])  # Move ax1 plot upwards



                    # Plot spectrogram
                    cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
                    ax2.set_xlabel('Time (s)')
                    f0lab = []
                    ax2.axvline(x=t0, c = '#377eb8', ls = '--', linewidth=0.7,label='Estimated arrival: '+str(np.round(t0,2))+' s')
                    
                    for pp in range(len(peaks)):
                        tprime = freqpeak[pp]
                        ft0p = peaks[pp]
                        f0 = calc_f0(tprime, t0, ft0p, v0, l, c)
                        
                        ft = calc_ft(times, t0, f0, v0, l, c)

                        ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) #(0,(5,10)),
                        
                        if np.abs(tprime -t0) < 1.5:
                            ax2.scatter(t0, ft0p, color='black', marker='x', s=30) 
                        f0lab.append(int(f0)) 
                        what_if = calc_ft(times, tarrive, fs/4, speed_mps, np.sqrt(dist_m**2 + alt_m**2), c)
                        #ax2.plot(times, what_if, 'red', ls = '--', linewidth=0.4)
                    f0lab_sorted = sorted(f0lab)
                    covm = np.sqrt(np.diag(covm))
                    if len(f0lab_sorted) <= 17:
                        fss = 'medium'
                    else:
                        fss = 'small'
                    ax2.set_title("Final Model:\nt0'= "+str(np.round(t0,2))+' +/- ' + str(np.round(covm[3],2)) + ' sec, v0 = '+str(np.round(v0,2))+' +/- ' + str(np.round(covm[1],2)) +' m/s, l = '+str(np.round(l,2))+' +/- ' + str(np.round(covm[2],2)) +' m, \n' + 'f0 = '+str(f0lab_sorted)+' +/- ' + str(np.round(covm[0],2)) +' Hz', fontsize=fss)
                    ax2.axvline(x=tarrive, c = '#e41a1c', ls = '--',linewidth=0.5,label='Wave arrvial: '+str(np.round(tarrive,2))+' s')

                    
                    ax2.legend(loc='upper right',fontsize = 'x-small')
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

                    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'
                    make_base_dir(BASE_DIR)
                    fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
                    plt.close()
                    
                    fig = plt.figure(figsize=(10,6))
                    plt.grid()

                    plt.plot(frequencies, arrive_time, c='#377eb8')
                    if auto_peak_pick == True:
                        peaks, _ = find_peaks(arrive_time, prominence=15)#0, distance = 10, height = 5, width=1) #for later change parameters for jets and permenant stations
                        np.diff(peaks)
                        plt.plot(peaks, arrive_time[peaks], c='k', marker="x", s=50)
                        for g in range(len(peaks)):
                            plt.text(peaks[g], arrive_time[peaks[g]], peaks[g], fontsize=15)
                    else:    
                        for pp in range(len(peaks)):
                            if np.abs(freqpeak[pp] -t0) < 1.5:
                                upper = int(peaks[pp] + 3)
                                lower = int(peaks[pp] - 3)
                                tt = spec[lower:upper, closest_index]
                                ampp = np.max(tt)
                                freqp = np.argmax(tt)+lower
                                plt.scatter(freqp, ampp, color='black', marker='x', s=100)
                                if isinstance(sta[n], int):
                                    plt.text(freqp - 5, ampp + 0.8, freqp, fontsize=17, fontweight='bold')
                                else:
                                    plt.text(freqp - 1, ampp + 0.8, freqp, fontsize=17, fontweight='bold')  
                    plt.xlim(0, int(fs/2))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    plt.ylim(0,vmax*1.1)
                    plt.xlabel('Frequency (Hz)', fontsize=17)
                    plt.ylabel('Relative Amplitude at t = {:.2f} s (dB)'.format(t0), fontsize=17)


                    make_base_dir('/scratch/irseppi/nodal_data/plane_info/5spec/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/')

                    fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(sta[n])+'_' + str(time[n]) + '.png')
                    plt.close() 
