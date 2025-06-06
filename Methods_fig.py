import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy.linalg as la
import obspy
import datetime
import pyproj
from scipy.signal import find_peaks, spectrogram
from prelude import *
from plot_func import *
import json

#Load seismometer data
seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

# Define the UTM projection with zone 6 and WGS84 ellipsoid
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

#Convert to UTM coordinates
seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]

seismo_utm_km = [(x, y) for x, y in zip(seismo_utm_x_km, seismo_utm_y_km)]

#Date, flight number and station used in this example
equip = "C185"
date = 20190214
flight = 528502194
sta = 1024

# Find find the data for the example station
for s in range(len(seismo_data)):
    if str(sta) == str(stations[s]):
        seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s]) 
        elevation = elevations[s] 
        break
    else:
        continue

#Load flight data
flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight) + '.csv'
flight_data = pd.read_csv(flight_file, sep=",")
flight_latitudes = flight_data['latitude']
flight_longitudes = flight_data['longitude']
time = flight_data['snapshot_id']
timestamp = flight_data['snapshot_id']
speed = flight_data['speed']
altitude = flight_data['altitude']

# Convert flight latitude and longitude to UTM coordinates
flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
flight_utm_x, flight_utm_y = zip(*flight_utm)

# Convert UTM coordinates to kilometers
flight_utm_x_km = [x / 1000 for x in flight_utm_x]
flight_utm_y_km = [y / 1000 for y in flight_utm_y]
flight_path = [(x,y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]

#Calling the function to find the closest point on the flight path to the seismometer
closest_p, dist_km, index = find_closest_point(flight_path, seismometer)
closest_x, closest_y = closest_p

# Calculate the approximate time of the closest approach
flight_utm_x1, flight_utm_y1 = flight_path[index]
flight_utm_x2, flight_utm_y2 = flight_path[index + 1]

# Calculate the difference vectors between flight coordinates
x_timestamp_dif_vec = flight_utm_x2 - flight_utm_x1
y_timestamp_dif_vec = flight_utm_y2 - flight_utm_y1

# Calculate the difference vectors between closest point and flight coordinates
cx_timestamp_dif_vec =  closest_x - flight_utm_x1
cy_timestamp_dif_vec = closest_y - flight_utm_y1

# Calculate the line vector and closest point vector magnitudes
line_vector = (x_timestamp_dif_vec, y_timestamp_dif_vec)
cline_vector = (cx_timestamp_dif_vec, cy_timestamp_dif_vec)

line_magnitude = np.sqrt(line_vector[0] ** 2 + line_vector[1] ** 2)
cline_magnitude = np.sqrt(cline_vector[0] ** 2 + cline_vector[1] ** 2)

# Calculate the length ratio and closest time
length_ratio = cline_magnitude / line_magnitude
closest_time = timestamp[index] + length_ratio*(timestamp[index+1] - timestamp[index])

# Convert altitude and speed to meters per second
alt_m = altitude[index] * 0.3048
speed_mps = speed[index] * 0.514444

# Calculate the height and distance in meters
height_m = alt_m - elevation 
dist_m = dist_km * 1000

# Set the midpoint time
tmid = closest_time

# Calculate the arrival time using the calc_time function
tarrive = calc_time(tmid,dist_m,height_m)

lon, lat = utm_proj(closest_x, closest_y, inverse=True)

input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'

file =  open(input_files, 'r')
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
            if abs(float(item['values'][i]) - float(alt_m)) < hold:
                hold = abs(float(item['values'][i]) - float(alt_m))
                z_index = i
folder_spec = equip + '_spec_c'
folder_spectrum = equip + '_spectrum_c'
for item in data_list:
    if item['parameter'] == 'T':
        Tc = - 273.15 + float(item['values'][z_index])

folder_spec =  equip + '_spec_cfc'
folder_spectrum = equip + '_spectrum_cfc'
c = speed_of_sound(Tc)
sound_speed = c

# Convert the arrival time to datetime object
ht = datetime.datetime.utcfromtimestamp(tarrive)
mins = ht.minute
secs = ht.second
h = ht.hour

h_u = str(h+1)
h = str(h)

# Read in seismic data
p = "/scratch/naalexeev/NODAL/2019-02-14T"+str(h)+":00:00.000000Z.2019-02-14T"+str(h_u)+":00:00.000000Z.1024.mseed"
tr = obspy.read(p)

# Trim the seismic 240 second time window
tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - 120, tr[2].stats.starttime + (mins * 60) + secs + 120)
data = tr[2][:]

# Create a title for the seismic data
title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'

# Get the time values of and sampling rate of the data
torg = tr[2].times()
fs = int(tr[2].stats.sampling_rate)

# Compute spectrogram
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
fig, ax = plt.subplots(6,1, figsize=(10, 8), sharex=True)
closest_index = np.argmin(np.abs(tarrive - times))
arrive_time = Sxx[:,closest_index]
for i in range(len(arrive_time)):
    if arrive_time[i] < 0:
        arrive_time[i] = 0

vmin = np.min(Sxx) 
vmax = np.max(Sxx)
ax[0].pcolormesh(times, frequencies, Sxx, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
# Calculate the median difference function (MDF)
a, b = Sxx.shape
MDF = np.zeros((a,b))
for row in range(len(Sxx)):
    m = len(Sxx[row])
    p = sorted(Sxx[row])
    median = p[int(m/2)]
    for col in range(m):
        MDF[row][col] = median

# Calculate the spectrogram with median normalization
spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

# Use the middle column of the spectrogram to intialize the minimum and maximum values for the color map
middle_index =  len(times) // 2
middle_column = spec[:, middle_index]
vmin = 0  
vmax = np.max(middle_column) 

# Set the velocity of the aircraft and distance between aircraft and station at the closest approach
tprime0 = 120
v0 = speed_mps
l = np.sqrt(dist_m**2 + (height_m)**2)
 

print("Please pick the points on the spectrogram that correspond to the primary overtone of the doppler curves.")

pick_again = 'y'
while pick_again == 'y':
    coords = []  # Reset the coordinates list
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

    # Function to handle mouse click events
    def onclick(event):
        global coords
        coords.append((event.xdata, event.ydata))  # Add clicked coordinates to the list
        plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add marker at clicked point
        plt.draw() 
        print('Clicked:', event.xdata, event.ydata)  

    # Connect the onclick function to the button press event
    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

    plt.show(block=True)  # Display the plot
    pick_again = input("Do you want to repick your points? (y or n)")  # Ask user if they want to repick

# Convert the list of coordinates to a numpy array
coords_array = np.array(coords)

# Initialize the initial model parameters
f0 = 116
m0 = [f0, v0, l, tprime0]

# Perform inversion using the initial model parameters and coordinates array
m,covm = invert_f(m0, coords_array, num_iterations=8)
f0 = m[0]
v0 = m[1]
l = m[2]
tprime0 = m[3]

# Calculate the theoretical arrival times
ft = calc_ft(times, tprime0, f0, v0, l, c)

# Find peaks in the middle column of the spectrogram
peaks = []
p, _ = find_peaks(middle_column, distance = 7)
corridor_width = (fs/2) / len(p)

# Adjust corridor width if no peaks are found
if len(p) == 0:
    corridor_width = fs/4

coord_inv = []

# Iterate over time values
for t_f in range(len(times)):
    upper = int(ft[t_f] + corridor_width)
    lower = int(ft[t_f] - corridor_width)
    if lower < 0:
        lower = 0
    if upper > len(frequencies):
        upper = len(frequencies)
    tt = spec[lower:upper, t_f]

    # Find the maximum amplitude frequency within the corridor
    max_amplitude_index = np.argmax(tt)
    max_amplitude_frequency = frequencies[max_amplitude_index+lower]
    peaks.append(max_amplitude_frequency)
    coord_inv.append((times[t_f], max_amplitude_frequency))

# Convert the list of coordinates to a numpy array
coord_inv_array = np.array(coord_inv)

# Perform inversion using the initial model parameters and coordinates array
m,_ = invert_f(m0, coord_inv_array, num_iterations=12)
f0 = m[0]
v0 = m[1]
l = m[2]
tprime0 = m[3]

# Calculate the theoretical arrival times
ft = calc_ft(times, tprime0, f0, v0, l, c)

# Calculate the difference between the theoretical arrival times and observed peaks
delf = np.array(ft) - np.array(peaks)

# Filter the coordinates based on the difference between theoretical and observed peaks
new_coord_inv_array = []
for i in range(len(delf)):
    if np.abs(delf[i]) <= 3:
        new_coord_inv_array.append(coord_inv_array[i])
coord_inv_array = np.array(new_coord_inv_array)

# Perform inversion using the filtered coordinates array
m,covm = invert_f(m0, coord_inv_array, num_iterations=12, sigma=5)
f0 = m[0]
v0 = m[1]
l = m[2]
tprime0 = m[3]

# Update the prior model parameters
mprior = []
mprior.append(v0)
mprior.append(l)
mprior.append(tprime0)

# Initialize the pick_again variable
pick_again = 'y'

print("Please pick one point on each overtone, it does not have to be at the center of the doppler.")

# Loop to allow user to repick points
while pick_again == 'y':
    peaks = []
    freqpeak = []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    plt.axvline(x=tprime0, c = '#377eb8', ls = '--')
    plt.axvline(x=120, c = '#e41a1c', ls = '--')

    # Function to handle mouse click events
    def onclick(event):
        global coords
        peaks.append(event.ydata)
        freqpeak.append(event.xdata)
        plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
        plt.draw() 
        print('Clicked:', event.xdata, event.ydata)  

    # Connect the onclick function to the button press event
    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

    plt.show(block=True)
    pick_again = input("Do you want to repick your points? (y or n)")


# Calculate the number of peaks
w = len(peaks)

# Iterate over each peak
for o in range(w):
    tprime = freqpeak[o]
    ft0p = peaks[o]
    
    # Calculate f0 using the given parameters
    f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
    
    # Append f0 to the prior model parameters
    mprior.append(f0)

# Convert mprior to a numpy array
mprior = np.array(mprior)

# Set the corridor width
corridor_width = 6

# Initialize lists for associated peaks, observed frequencies, and observed times
peaks_assos = []
fobs = []
tobs = []

# Iterate over each peak
for pp in range(len(peaks)):
    tprime = freqpeak[pp]
    ft0p = peaks[pp]
    
    # Calculate f0 using the given parameters
    f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
    
    # Initialize lists for maximum frequencies, coordinates, and times
    maxfreq = []
    coord_inv = []
    ttt = []

    # Calculate the upper and lower frequencies within the corridor
    f01 = f0 + corridor_width
    f02 = f0 - corridor_width
    upper = calc_ft(times, tprime0, f01, v0, l, c)
    lower = calc_ft(times, tprime0, f02, v0, l, c)

    # Iterate over each time
    for t_f in range(len(times)):
        try:
            # Extract the corresponding spectrogram values within the frequency corridor
            tt = spec[int(np.round(lower[t_f], 0)):int(np.round(upper[t_f], 0)), t_f]
            
            try:
                # Find the maximum amplitude frequency within the corridor
                max_amplitude_index, _ = find_peaks(tt, prominence=15, wlen=10, height=vmax*0.1)
                maxa = np.argmax(tt[max_amplitude_index])
                max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa]) + int(np.round(lower[t_f], 0))]
            except:
                continue
            
            # Append the maximum amplitude frequency and its coordinates
            maxfreq.append(max_amplitude_frequency)
            coord_inv.append((times[t_f], max_amplitude_frequency))
            ttt.append(times[t_f])

        except:
            continue
    
    if len(coord_inv) > 0:
        if f0 < 200:
            # Perform inversion using the filtered coordinates array and updated model parameters
            coord_inv_array = np.array(coord_inv)
            mtest = [f0, v0, l, tprime0]
            mtest, _ = invert_f(mtest, coord_inv_array, num_iterations=4)
            ft = calc_ft(ttt, mtest[3], mtest[0], mtest[1], mtest[2], c)
        else:
            # Calculate the theoretical arrival times using the original model parameters
            ft = calc_ft(ttt, tprime0, f0, v0, l, c)

        # Calculate the difference between the theoretical arrival times and observed maximum frequencies
        delf = np.array(ft) - np.array(maxfreq)

        count = 0
        # Iterate over each difference
        for i in range(len(delf)):
            if np.abs(delf[i]) <= 3:
                # Append the associated peak and observed frequency
                fobs.append(maxfreq[i])
                tobs.append(ttt[i])
                count += 1
        
        # Append the count of associated peaks
        peaks_assos.append(count)
    else:
        continue

print('Please pick two points on the spectrogram that correspond to the start and end of the time window you want pull data from in the inversion.')

# Initialize variable to prompt for repicking points
pick_again = 'y'

# Loop to allow for repicking points
while pick_again == 'y':
    set_time = []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    plt.scatter(tobs,fobs, color='black', marker='x')
    
    # Function to handle mouse click event
    def onclick(event):
        global coords
        set_time.append(event.xdata) 
        plt.scatter(event.xdata, event.ydata, color='red', marker='x')  # Add this line
        plt.draw() 
        print('Clicked:', event.xdata, event.ydata)  
    
    # Connect the mouse click event to the function
    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show(block=True)

    # Prompt user to repick points
    pick_again = input("Do you want to repick your points? (y or n)")

# Set the start and end time based on the picked points
start_time = set_time[0]
end_time = set_time[1]
ftobs = []
ffobs = []

ftobs = []

peak_ass = []
cum = 0
for p in range(w):
    count = 0
    for j in range(cum,cum+peaks_assos[p]):
        if tobs[j] >= start_time and tobs[j] <= end_time:
            ftobs.append(tobs[j])
            ffobs.append(fobs[j])
            count += 1
    cum = cum + peaks_assos[p]

    peak_ass.append(count)
peaks_assos = peak_ass
tobs = ftobs
fobs = ffobs

# Initialize variables
qv = 0
num_iterations = 4

# Define prior covariance matrix
cprior = np.zeros((w+3,w+3))
for row in range(len(cprior)):
    if row == 0:
        cprior[row][row] = 20**2
    elif row == 1:
        cprior[row][row] = 500**2
    elif row == 2:
        cprior[row][row] = 20**2
    else:
        cprior[row][row] = 1**2

# Define data covariance matrix
Cd = np.zeros((len(fobs), len(fobs)), int)
np.fill_diagonal(Cd, 3**2)

# Initialize model parameters
mnew = np.array(mprior)

# Perform iterations
while qv < num_iterations:
    G = np.zeros((0,w+3))
    fnew = []
    cum = 0
    for p in range(w):
        new_row = np.zeros(w+3)
        tprime = freqpeak[p]
        ft0p = peaks[p]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)

        # Calculate derivatives and update G matrix
        for j in range(cum,cum+peaks_assos[p]):
            tprime = tobs[j]
            t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
            ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

            f_derivef0, f_derivev0, f_derivel, f_derivetprime0 = df(f0,v0,l,tprime0, tobs[j])
        
            new_row[0] = f_derivev0
            new_row[1] = f_derivel
            new_row[2] = f_derivetprime0
            new_row[3+p] = f_derivef0
                    
            G = np.vstack((G, new_row))
                    
            fnew.append(ft0p)
    
        cum = cum + peaks_assos[p]

    # Update model parameters using least squares inversion
    m = np.array(mnew) + cprior@G.T@la.inv(G@cprior@G.T+Cd)@(np.array(fobs)- np.array(fnew))
    mnew = m
    v0 = mnew[0]
    l = mnew[1]
    tprime0 = mnew[2]
    f0_array = mnew[3:]

    print(m)
    qv += 1

# Calculate covariance matrix
covm = la.inv(G.T@la.inv(Cd)@G + la.inv(cprior))

# Find the index of the closest time value to tprime0
closest_index = np.argmin(np.abs(tprime0 - times))

# Get the arrival time values at the closest time index
arrive_time = spec[:, closest_index]

# Set negative arrival time values to 0
for i in range(len(arrive_time)):
    if arrive_time[i] < 0:
        arrive_time[i] = 0

# Calculate the minimum and maximum values of the arrival time
vmin = np.min(arrive_time)
vmax = np.max(arrive_time)


# Create subplots for the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

# Plot the original data
ax1.plot(torg, data, 'k', linewidth=0.5)
ax1.set_title(title)

ax1.margins(x=0)
ax1.set_position([0.125, 0.6, 0.775, 0.3])  # Move ax1 plot upwards

# Plot spectrogram
cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
ax2.set_xlabel('Time (s)')

# Plot estimated arrival time
ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label='Estimated arrival: '+str(np.round(tprime0,2))+' s')

# Calculate and plot the model
covm = np.sqrt(np.diag(covm))
for pp in range(len(f0_array)):
    f0 = f0_array[pp]
    
    ft = calc_ft(times, tprime0, f0, v0, l, c)

    ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 
    tprime = tprime0
    t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
    ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
    
    ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 

# Set the title for the final model plot
fss = 'x-small'
f0lab = sorted(f0_array)
for i in range(len(f0lab)):
    f0lab[i] = (str(np.round(f0lab[i],2)))
ax2.set_title("Final Model:\nt0'= "+str(np.round(tprime0,2)) + ' sec, v0 = '+str(np.round(v0,2)) +' m/s, l = '+str(np.round(l,2)) +' m, \n' + 'f0 = '+', '.join(f0lab) +' Hz', fontsize=fss)

# Plot wave arrival line
ax2.axvline(x=120, c = '#e41a1c', ls = '--',linewidth=0.5,label='Wave arrvial: 120 s')

# Add legend and labels to ax2
ax2.legend(loc='upper right',fontsize = 'x-small')
ax2.set_ylabel('Frequency (Hz)')

ax2.margins(x=0)
ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

# Add colorbar
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

# Display the figure
plt.show()     
