import pyproj
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from pathlib import Path
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks
from prelude import make_base_dir, calc_ft, calc_f0, invert_f

################################################################################################################################################

def plot_map(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km, closest_time, d, index, flight_num, date, seismometer, closest_p, head, station):
    """
    Plot the flight path and seismometer locations on a map.

    Args:
        flight_utm_x_km (list): List of UTM x coordinates for the flight path.
        flight_utm_y_km (list): List of UTM y coordinates for the flight path.
        seismo_utm_x_km (list): List of UTM x coordinates for the seismometers.
        seismo_utm_y_km (list): List of UTM y coordinates for the seismometers.
        closest_time (str): The closest time to the seismometer.
        d (float): Distance to the closest point in kilometers.
        index (int): Index of the closest point in the flight path.
        flight_num (int): Flight number.
        date (str): Date of the flight.
        seismometer (tuple): Coordinates of the seismometer.
        closest_p (tuple): Coordinates of the closest point on the flight path.
        head (list): List of headings for the flight path.
        station (str): Station identifier.

    Returns:
        None
    """

    utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

    closest_x, closest_y = closest_p
    min_lon = -150.7
    max_lon = -147.3
    min_lat = 62.2
    max_lat = 65.3
    lxmin,lymin = utm_proj(min_lon, min_lat)
    lxmax, lymax = utm_proj(max_lon, max_lat)
    y = [closest_y, seismometer[1]]
    x = [closest_x, seismometer[0]]
    yy = sum(y) / len(y)
    xx = sum(x) / len(x)
    
    #set the size of the map depending on the distance
    min_x = int(xx - 2)
    max_x = int(xx + 2)
    min_y = int(yy - 2)
    max_y = int(yy + 2)
    if d < 0.5:
        min_x = (xx - 1)
        max_x = (xx + 1)
        min_y = (yy - 1)
        max_y = (yy + 1)
    if d < 0.1:
        min_x = (xx - 0.1)
        max_x = (xx + 0.1)
        min_y = (yy - 0.1)
        max_y = (yy + 0.1)

    # Create a figure with two subplots side by side
    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 2]}) 
    fig.subplots_adjust(wspace=1) 

    # Ploting the entire array and the flight path
    axs[0].set_xticks(np.arange(int(lxmin/1000)-7, int(lxmax/1000), 50))
    axs[0].set_yticks(np.arange(int(lymin/1000)-1, int(lymax/1000), 50))
    axs[0].set_xlabel('UTM Easting (km)')
    axs[0].set_ylabel('UTM Northing (km)')
    axs[0].set_aspect('equal')
    axs[1].set_aspect('equal')

    axs[0].grid(True, linestyle='dotted', color='gray')
    axs[1].grid(True, linestyle='dotted', color='gray')

    axs[0].scatter(seismo_utm_x_km, seismo_utm_y_km, c='#e41a1c', s = 3, label='seismometers')
    axs[0].plot(flight_utm_x_km, flight_utm_y_km, '-', c='#377eb8', lw=1, ms = 1, label='flight path')

    # Plot arrows indicating the direction of the flight path
    for i in range(1, len(flight_utm_y_km)-1, int(len(flight_utm_y_km)/5)):
        direction = np.arctan2(flight_utm_y_km[i+1] - flight_utm_y_km[i], flight_utm_x_km[i+1] - flight_utm_x_km[i])
        if (flight_utm_x_km[i+1] - flight_utm_x_km[i]) == 0:
            continue
        m = (flight_utm_y_km[i+1] - flight_utm_y_km[i])/(flight_utm_x_km[i+1] - flight_utm_x_km[i])
        if m == 0:
            continue
        b = flight_utm_y_km[i] - m*flight_utm_x_km[i]
        axs[0].quiver((flight_utm_y_km[i]-b)/m, flight_utm_y_km[i], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', headwidth = 10)

    # Set labels and title
    axs[0].set_xlim(int(lxmin/1000), int(lxmax/1000))
    axs[0].set_ylim(int(lymin/1000), int(lymax/1000))
    axs[0].tick_params(axis='both', which='major', labelsize=9)

    head_avg = (head[index]+head[index+1])/2
    converted_angle = (90 - head_avg) % 360
    heading = np.deg2rad(converted_angle)

    # Define the UTM and latitude/longitude coordinate systems
    rect = Rectangle((min_x, min_y), (max_x-min_x), (max_y-min_y), ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
    axs[0].add_patch(rect)

    axs[1].set_xticks(np.arange(min_x, max_x, np.round(((max_x - min_x) / 4), 1)))
    axs[1].set_yticks(np.arange(min_y, max_y, np.round(((max_y - min_y) / 4), 1)))

    # Plot the zoomed in map on the second subplot
    axs[1].plot(x,y, '--', c='#ff7f00')
    axs[1].plot(flight_utm_x_km, flight_utm_y_km, c='#377eb8',linestyle ='dotted')
    axs[1].scatter(flight_utm_x_km, flight_utm_y_km, c='#377eb8', s=20)

    axs[1].set_xlim(min_x, max_x)
    axs[1].set_ylim(min_y, max_y)

    axs[1].tick_params(axis='both', which='major', labelsize=9)
    axs[1].ticklabel_format(useOffset=False, style='plain')

    direction = np.arctan2(flight_utm_y_km[index+1] - flight_utm_y_km[index], flight_utm_x_km[index+1] - flight_utm_x_km[index])
    
    axs[1].quiver(closest_x, closest_y, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=8)

    axs[1].quiver(closest_x, closest_y, np.cos(heading), np.sin(heading), angles='xy', scale = 8, color='#999999')
    axs[1].scatter(closest_x, closest_y, c='#377eb8', s=50, zorder=3)
    axs[1].scatter(seismometer[0], seismometer[1], c='#e41a1c', s=50, zorder=3)

    axs[1].text(xx,yy, str(round(d, 2))+' km', fontsize=12, fontweight='bold')

    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
    con = mpatch.ConnectionPatch(xyA=(max_x, min_y), xyB=(min_x, min_y), coordsA="data", coordsB="data", axesA=axs[0], axesB=axs[1], color="black", linestyle="--")
    fig.add_artist(con)
    con = mpatch.ConnectionPatch(xyA=(max_x, max_y), xyB=(min_x, max_y), coordsA="data", coordsB="data", axesA=axs[0], axesB=axs[1], color="black", linestyle="--")
    fig.add_artist(con)
    plt.tight_layout(pad=0.2, rect=[0, 0, 1, 1])

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/map_all_UTM/' + str(date) + '/'+ str(flight_num) + '/' + str(station) + '/'
    make_base_dir(BASE_DIR)
    plt.savefig('/scratch/irseppi/nodal_data/plane_info/map_all_UTM/' + str(date) + '/' + str(flight_num) + '/' + str(station) + '/map_' + str(flight_num) + '_' + str(closest_time) + '.png', bbox_inches='tight')
    plt.close()


##############################################################################################################################################################################################################

def remove_median(Sxx):
    """
    Remove the median from the spectrogram.

    Args:
        Sxx (array): The spectrogram data.

    Returns:
        spec: The spectrogram data with the median removed
        MDF: The median removed from the spectrogram
    """

    a, b = Sxx.shape

    MDF = np.zeros((a,b))
    for row in range(len(Sxx)):
        median = np.median(Sxx[row])
        MDF[row, :] = median

    # Avoid log10(0) by replacing zeros with a small positive value
    Sxx_safe = np.where(Sxx == 0, 1e-10, Sxx)
    MDF_safe = np.where(MDF == 0, 1e-10, MDF)

    spec = 10 * np.log10(Sxx_safe) - (10 * np.log10(MDF_safe))
    return spec, MDF

############################################################################################################################################################################################################################

def plot_spectrogram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, Cpost, flight, middle_index, tarrive, closest_time, dir_name, plot_show=True):
    """
    Plot and save the waveform, unfiltered, and the spectrogram of the given data. Include the estimated curve using the final model parameters outputs from the inversions and tprime0 initial guess compared to the final.

    Args:
        data (array): The waveform data.
        fs (int): The sampling frequency.
        torg (array): The time array.
        title (str): The title of the plot.
        spec (array): The spectrogram data.
        times (array): The time array for the spectrogram.
        frequencies (array): The frequency array for the spectrogram.
        tprime0 (float): The estimated arrival time.
        v0 (float): The velocity.
        l (float): The distance.
        c (float): The speed of sound.
        f0_array (array): The array of frequencies.
        F_m (float): The data misfit value.
        arrive_time (array): The arrival time array.
        MDF (array): Median removed from spectrogram.
        Cpost (array): The normalized posterior covariance matrix.
        flight (int): The flight number.
        middle_index (int): The index of the middle column.
        tarrive (float): The initial estimation of the arrival time at the station from wave generated at aircraft at closest approach.
        closest_time (float): The time of closest approach of aircraft, for saving the file.
        dir_name (str): The directory name.
        plot_show (bool): If True, show the plot and ask user to provide a quality number. If False, save the plot without showing it. 

    Returns:
        str: The user assigned quality number.
    """
    # Plot settings and calculations
    vmin = np.min(arrive_time) 
    vmax = np.max(arrive_time)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

    ax1.plot(torg, data, 'k', linewidth=0.5)
    ax1.set_title(title)

    ax1.margins(x=0)
    ax1.set_position([0.125, 0.6, 0.775, 0.3])  # Move ax1 plot upwards
    ax1.axvline(x=tarrive, c = '#e41a1c', ls = '--',linewidth=0.5,label= r'$t_{i}$ = ' + "%.2f" % tarrive +' s')
    ax1.axvline(x=tprime0, c = 'k', ls = '--',linewidth=0.5)
    # Plot spectrogram
    cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
    ax2.set_xlabel('Time (s)')
    f0lab = []
    ax2.axvline(x=tarrive, c = '#e41a1c', ls = '--',linewidth=0.5,label= r'$t_{i}$ = ' + "%.2f" % tarrive +' s')

    ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label= "t\u2080' = " + "%.2f" % tprime0 +' s')
    for pp in range(len(f0_array)):
        f0 = f0_array[pp]
        
        ft = calc_ft(times, tprime0, f0, v0, l, c)

        ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 
        tprime = tprime0
        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
        
        ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 

    fss = 'x-small'
    f0lab = sorted(f0_array)

    if len(f0_array) <= 1:
        med_df = "NaN"
        mad_df = "NaN"

    else:
        #Generate random samples of f0 values withing their sigma from the covariance matrix 
        #Calculate the median of the differences and MAD to obtain error
        f_range = []
        NTRY = 1000
        for N in range(NTRY):
            ftry = []
            for c_index  in range(4, len(Cpost)):
                xmin = f0_array[c_index-4] - Cpost[c_index]
                xmax = f0_array[c_index-4] + Cpost[c_index]
                xtry = xmin + (xmax-xmin)*np.random.rand()
                ftry.append(xtry)

            ftry = np.sort(ftry)
            f1 = []
            for g in range(len(ftry)):
                if g == 0:
                    continue
                diff = ftry[g] - ftry[g - 1]
                f1.append(diff)
            med = np.nanmedian(f1)
            f_range.append(med)
        med_df = np.nanmedian(f_range)
        mad_df = np.nanmedian(np.abs(f_range - med_df))
    if med_df == "NaN":
        ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % Cpost[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % Cpost[0]+' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % Cpost[1] + ' m, \n' + 'f\u2080 = ['+', '.join(["%.2f" % f for f in f0lab]) +'] \u00B1 ' + "%.2f" % np.median(Cpost[3:]) +' Hz, df\u2080 = NaN \u00B1 NaN Hz\nMisfit: ' + "%.4f" % F_m, fontsize=fss)
    else:
        ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % Cpost[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % Cpost[0] +' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % Cpost[1] + ' m, \n' + 'f\u2080 = ['+', '.join(["%.2f" % f for f in f0lab]) +'] \u00B1 ' + "%.2f" % np.median(Cpost[3:]) +' Hz, df\u2080 = ' + "%.2f" % med_df + ' \u00B1 ' + "%.2f" % mad_df + ' Hz\nMisfit: ' + "%.4f" % F_m + ' [FH/VT]', fontsize=fss)

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
    if plot_show:
        plt.show()     
        qnum = input('What quality number would you give this?(first num for data quality(0-3), second for ability to fit model to data(0-1))')
    else:
        qnum = '__'
    
    fig.savefig(dir_name+'/'+str(closest_time)+'_'+str(flight)+'.png')
    plt.close()

    return qnum

################################################################################################################################################################################################################################################################################################################################

def plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, dir_name):
    """
    Plot and save the spectrum with markers for the overtones arriving at the station at tprime0.

    Args:
        spec (numpy.ndarray): The spectrum data.
        frequencies (numpy.ndarray): The frequencies.
        tprime0 (float): The reference arrival time.
        v0 (float): The velocity.
        l (float): The distance.
        c (float): The speed of light.
        f0_array (list): The list of frequencies.
        arrive_time (numpy.ndarray): The arrival times.
        fs (int): The sampling frequency.
        closest_index (int): The index of the closest time.
        closest_time (float): The closest time.
        sta (int or str): The station identifier.
        dir_name (str): The directory name.

    Returns:
        None
    """
    vmax = np.max(arrive_time)
    fig = plt.figure(figsize=(10,6))
    plt.grid()

    plt.plot(frequencies, spec[:,closest_index], c='#377eb8')
        
    for pp in range(len(f0_array)):
        f0 = f0_array[pp]
        if fs/2 < f0:
            continue
        tprime = tprime0
        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
        if ft0p == np.nan:
            continue
        if ft0p > 250:
            continue
        
        upper = int(ft0p + 10)
        lower = int(ft0p - 10)
        tt = spec[lower:upper, closest_index]
        if upper > 250:
            freqp = ft0p
            ampp = np.interp(ft0p, frequencies, arrive_time)
        elif lower < 0:
            freqp = ft0p
            ampp = np.interp(ft0p, frequencies, arrive_time)
        else:
            ampp = np.max(tt)
            freqp = np.argmax(tt)+lower
        plt.scatter(freqp, ampp, color='black', marker='x', s=100)
        if isinstance(sta, int):
            plt.text(freqp - 5, ampp + 0.8, freqp, fontsize=17, fontweight='bold')
        else:
            plt.text(freqp - 1, ampp + 0.8, freqp, fontsize=17, fontweight='bold')  

    plt.xlim(0, int(fs/2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0,vmax*1.1)
    plt.xlabel('Frequency (Hz)', fontsize=17)
    plt.ylabel('Relative Amplitude at t = {:.2f} s (dB)'.format(tprime0), fontsize=17)

    fig.savefig(dir_name + '/'+str(sta)+'_' + str(closest_time) + '.png')
    plt.close()

##############################################################################################################################################################################################################

def doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight, sta, equip, closest_time, tarrive, make_picks=True, spec_window = 120):
    """
    Pick the points for the doppler shift.

    Args:
        spec (numpy.ndarray): The spectrogram data.
        times (numpy.ndarray): The time array.
        frequencies (numpy.ndarray): The frequency array.
        vmin (float): The minimum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        vmax (float): The maximum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        month (int): The month of the data.
        day (int): The day of the data.
        flight (int): The flight number.
        sta (int or str): The station identifier.
        equip (str): The equipment identifier.
        closest_time (float): The time of closest approach.
        start_time (float): The start time of the spectrogram, to save for future reference on plotting the spectrogram.
        make_picks (bool): If you come to this function and no picks exist, it will allow you to make new picks.

    Returns:
        list: The list of picks the user picked along the most prominent overtone.
    """

    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/inversepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight) + '.csv'

    if Path(file_name).exists():
        coords = []
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                coords.append((float(pick_data[0]), float(pick_data[1])))
            if len(pick_data) == 4:
                start_time = float(pick_data[2])
            else:
                plt.figure()
                plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                plt.scatter([coord[0] for coord in coords], [coord[1] for coord in coords], color='black', marker='x')
                plt.show()
                correct_time = input("No start time found in file. Do your picks line up with this signal?(y/n): ")
                if correct_time == 'y': 
                    start_time = tarrive - spec_window
                    # Rewrite file with start_time as third column and move \n to next column
                    with open(file_name, 'r') as file:
                        lines = file.readlines()
                    with open(file_name, 'w') as file:
                        for line in lines:
                            pick_data = line.strip().split(',')
                            # Only keep first two columns, append start_time, then move \n to next column
                            if len(pick_data) >= 2:
                                file.write(f"{pick_data[0]},{pick_data[1]},{start_time},\n")
                else:
                    return [], None
        file.close()  
        return coords, start_time
    
    elif make_picks:
        BASE_DIR = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/inversepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight) + '/' + str(sta) + '/'
        make_base_dir(BASE_DIR)
        pick_again = 'y'
        start_time = tarrive - spec_window
        while pick_again == 'y':
            r1 = open(file_name, 'w')
            coords = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            def onclick(event, coords=coords):
                #global coords
                coords.append((event.xdata, event.ydata))
                plt.scatter(event.xdata, event.ydata, color='black', marker='x')  
                plt.draw() 
                print('Clicked:', event.xdata, event.ydata)  
                r1.write(str(event.xdata) + ',' + str(event.ydata) + ',' + str(start_time) + ',\n')
            plt.gcf().canvas.mpl_connect('button_press_event', onclick)

            plt.show(block=True)
            r1.close()
            pick_again = input("Do you want to repick your points? (y or n)")
        return coords, start_time
    else:
        return [], None

##############################################################################################################################################################################################################

def overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight, sta, equip, closest_time, start_time, tprime0, tarrive, make_picks=True):
    """
    Pick the points for the overtone shift.

    Args:
        spec (numpy.ndarray): The spectrogram data.
        times (numpy.ndarray): The time array.
        frequencies (numpy.ndarray): The frequency array.
        vmin (float): The minimum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        vmax (float): The maximum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        month (int): The month of the data.
        day (int): The day of the data.
        flight (int): The flight number.
        sta (int or str): The station identifier.
        equip (str): The equipment identifier.
        closest_time (float): The time of closest approach.
        start_time (float): The start time of the spectrogram, to save for future reference on plotting the spectrogram.
        tprime0 (float): The estimated acoustic wave arrival time.
        tarrive (float): The initial calculated time of acoustic wave arrival.
        make_picks (bool): If you come to this function and no picks exist, it will allow you to make new picks.

    Returns:
        list: List of frequencies picked by user along different overtones.
        list: List of times corresponding to the picked frequencies.
    """

    output2 = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/overtonepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight) + '.csv'
    if Path(output2).exists():

        peaks = []
        freqpeak = []
        with open(output2, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                peaks.append(float(pick_data[1]))
                freqpeak.append(float(pick_data[0]))
        file.close()  
        return peaks, freqpeak
    
    elif make_picks:
        BASE_DIR = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/overtonepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight) + '/' + str(sta) + '/'
        make_base_dir(BASE_DIR)
        pick_again = 'y'
        while pick_again == 'y':
            r2 = open(output2, 'w')
            peaks = []
            freqpeak = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            plt.axvline(x=tprime0, c = '#377eb8', ls = '--')
            plt.axvline(x=tarrive-start_time, c = '#e41a1c', ls = '--')
            def onclick(event):
                #global coords
                peaks.append(event.ydata)
                freqpeak.append(event.xdata)
                plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                plt.draw() 
                print('Clicked:', event.xdata, event.ydata)  
                r2.write(str(event.xdata) + ',' + str(event.ydata) + ',' + str(start_time) + ',\n')
            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

            plt.show(block=True)
            r2.close()
            pick_again = input("Do you want to repick you points? (y or n)")
        
        return peaks, freqpeak
    else:
        return [], []

##############################################################################################################################################################################################################

def time_picks(month, day, flight, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, w, peaks_assos, make_picks=True):
    """
    Pick the points for the time shift.

    Args:

        month (int): The month of the data.
        day (int): The day of the data.
        flight (int): The flight number.
        sta (int or str): The station identifier.
        equip (str): The equipment identifier.
        tobs (list): The time array.
        fobs (list): The frequency array.
        closest_time (float): The time of closest approach.
        start_time (float): The start time of the spectrogram, to save for future refrence on plotting the spectrogram.
        spec (numpy.ndarray): The spectrogram data.
        times (numpy.ndarray): The time array.
        frequencies (numpy.ndarray): The frequency array.
        vmin (float): The minimum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        vmax (float): The maximum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        w (int): The number of peaks.
        peaks_assos (list or bool): The number of peaks associated with each overtone.
        make_picks (bool): If you come to this function and no picks exist, it will allow you to make new picks.

    Returns:
        list: The time array, including data for all overtones.
        list: The frequency array, including data for all overtones.
        list: The number of data points associated with each overtone, for indexing purposes.
    """

    output3 = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/timepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight) + '.csv'
    if Path(output3).exists():
        set_time = []
        with open(output3, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                set_time.append(float(pick_data[0]))
        file.close()  
        if len(set_time) == 0:
            return tobs, fobs, peaks_assos
        s_time = set_time[0]
        e_time = set_time[1]
        ftobs = []
        ffobs = []
     
        peak_ass = []
        cum = 0
        
        for p in range(w):
            count = 0
            for j in range(cum,cum+peaks_assos[p]):
                if tobs[j] >= s_time and tobs[j] <= e_time:
                    ftobs.append(tobs[j])
                    ffobs.append(fobs[j])
                    count += 1
            cum = cum + peaks_assos[p]
        
            peak_ass.append(count)
        peaks_assos = peak_ass
        tobs = ftobs
        fobs = ffobs

        return tobs, fobs, peaks_assos

    elif make_picks:
        BASE_DIR = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/timepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'
        make_base_dir(BASE_DIR)
        
        pick_again = 'y'
        while pick_again == 'y':
            r3 = open(output3,'w')
            set_time = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            plt.scatter(tobs,fobs, color='black', marker='x')
            def onclick(event):
                #global coords
                set_time.append(event.xdata) 
                plt.scatter(event.xdata, event.ydata, color='red', marker='x')  # Add this line
                plt.draw() 
                print('Clicked:', event.xdata, event.ydata)  
                r3.write(str(event.xdata) + ',' + str(event.ydata) + ',' + str(start_time) + ',\n')

            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
            plt.show(block=True)
            plt.close()
            r3.close()
            pick_again = input("Do you want to repick you points? (y or n)")
        if len(set_time) == 0:
            return tobs, fobs, peaks_assos
        s_time = set_time[0]
        e_time = set_time[1]
        ftobs = []
        ffobs = []

        peak_ass = []
        cum = 0
        for p in range(w):
            count = 0
            for j in range(cum,cum+peaks_assos[p]):
                if tobs[j] >= s_time and tobs[j] <= e_time:
                    ftobs.append(tobs[j])
                    ffobs.append(fobs[j])
                    count += 1
            cum = cum + peaks_assos[p]
        
            peak_ass.append(count)
        peaks_assos = peak_ass
        tobs = ftobs
        fobs = ffobs

        return tobs, fobs, peaks_assos
    else:
        return tobs, fobs, peaks_assos

###############################################################################################################################

def get_auto_picks_1o(times, frequencies, spec, ft, corridor_width, mprior, sigma_prior):
    """
    Get automatic picks for the first overtone.

    Args:
        times (np.ndarray): Array of time values.
        frequencies (np.ndarray): Array of frequency values.
        spec (np.ndarray): Spectrogram data.
        ft (np.ndarray): Frequency calculated from model parameters using `calc_ft`.
        corridor_width (float): Width of the corridor for picking.

    Returns:
        np.ndarray: Array of auto picked coordinates.

    """
    coord_inv = []

    for t_f in range(len(times)):

        upper = int(ft[t_f] + corridor_width)
        lower = int(ft[t_f] - corridor_width)

        if lower < 0:
            lower = 0
        elif lower >= 250:
            continue
        else:
            pass
        if upper > 250:
            upper = 250

        tt = spec[lower:upper, t_f]
        max_amplitude_index = np.argmax(tt)

        max_amplitude_frequency = frequencies[max_amplitude_index+lower]
        coord_inv.append((times[t_f], max_amplitude_frequency))

    coord_inv_array = np.array(coord_inv)

    m,_,_,_ = invert_f(mprior,sigma_prior, coord_inv_array,num_iterations=5)
    f0 = m[0]
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]

    ft = calc_ft(coord_inv_array[:, 0], tprime0, f0, v0, l, c)

    delf = np.array(ft) - np.array(coord_inv_array[:, 1])
    
    new_coord_inv_array = []
    for i in range(len(delf)):
        if np.abs(delf[i]) <= 3:
            new_coord_inv_array.append(coord_inv_array[i])
    coord_inv_array = np.array(new_coord_inv_array)

    return coord_inv_array

################################################################################################################################

def get_auto_picks_full(peaks, time_peaks, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax, equip):
    """
    Get automatic picks for all overtones.

    Args:
        peaks (list): List of peak frequencies.
        time_peaks (list): List of times corresponding to the peaks.
        times (np.ndarray): Array of time values from fft.
        frequencies (np.ndarray): Array of frequency values from fft.
        spec (np.ndarray): Spectrogram data from fft.
        corridor_width (float): Width of the corridor for picking.
        tprime0 (float): Model parameter for the arrival time.
        v0 (float): Model parameter for the velocity.
        l (float): Model parameter for the distance.
        c (float): Model parameter for the speed of sound.
        sigma_prior (float): Prior uncertainty for the model parameters.
        vmax (float): Maximum amplitude value for peak detection.
        equip (str): Equipment identifier.

    Returns:
        list: List of observed times.
        list: List of observed frequencies.
        list: List of counts of peaks associated with each overtone, for indexing.
        list: List of fundamental frequencies calculated for each peak.
    """

    peaks_assos = []
    fobs = []
    tobs = []
    f0_array = []
  
    for pp in range(len(peaks)):
        tprime = time_peaks[pp]
        ft0p = peaks[pp]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        f0_array.append(f0)

        maxfreq = []
        coord_inv = []
        ttt = []

        ft = calc_ft(times,  tprime0, f0, v0, l, c)

        for t_f in range(len(times)):

            upper = int(ft[t_f] + corridor_width)
            lower = int(ft[t_f] - corridor_width)

            if lower < 0:
                lower = 0
            elif lower >= 250:
                continue
            else:
                pass
            if upper > 250:
                upper = 250
            
            tt = spec[int(np.round(lower,0)):int(np.round(upper,0)), t_f]

            max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
            if len(max_amplitude_index) == 0:
                continue

            maxa = np.argmax(tt[max_amplitude_index])
            max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa])+int(np.round(lower,0))]

            maxfreq.append(max_amplitude_frequency)
            coord_inv.append((times[t_f], max_amplitude_frequency))
            ttt.append(times[t_f])

        if len(ttt) > 0 and f0 <= 230:
            coord_inv_array = np.array(coord_inv)
            mtest = [f0,v0, l, tprime0,c]
            mtest,_,_,_ = invert_f(mtest,sigma_prior, coord_inv_array, num_iterations=4)
            ft = calc_ft(ttt,  mtest[3], mtest[0], mtest[1], mtest[2], mtest[4])
            delf = np.array(ft) - np.array(maxfreq)

            count = 0
            for i in range(len(delf)):
                if np.abs(delf[i]) <= (4):
                    fobs.append(maxfreq[i])
                    tobs.append(ttt[i])
                    count += 1
            peaks_assos.append(count)
        elif f0 > 230:
            for i in range(len(ttt)):
                fobs.append(maxfreq[i])
                tobs.append(ttt[i])
            peaks_assos.append(len(maxfreq))
        else:
            peaks_assos.append(0)

    return tobs, fobs, peaks_assos, f0_array
