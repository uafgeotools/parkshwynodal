import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy.linalg as la
import obspy
import datetime
import pyproj
from scipy.signal import find_peaks, spectrogram
from src.main_inv_fig_functions import remove_median, get_auto_picks_full, get_auto_picks_1o
from src.doppler_funcs import invert_f, calc_ft, full_inversion, calc_f0
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from matplotlib.ticker import MaxNLocator
# Initialize a client for a specific FDSN data center (e.g., IRIS, GEONET)
client = Client("http://service.iris.edu", service_mappings={"dataselect": "http://service.iris.edu/ph5ws/dataselect/1"}) 

starttime = UTCDateTime("2019-03-04T01:17:22")
endtime = UTCDateTime("2019-03-04T01:21:22")

st = client.get_waveforms("ZE", "1010", "*", "DPZ", starttime, endtime)
tr = st[0]

data = tr.data

# Create a title for the seismic data
title = f'{tr.stats.network}.{tr.stats.station}.{tr.stats.location}.{tr.stats.channel} − starting {tr.stats["starttime"]}'

# Get the time values of and sampling rate of the data
torg = tr.times()
fs = int(tr.stats.sampling_rate)

# Compute spectrogram
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 

# Calculate the median difference function (MDF)
spec, MDF = remove_median(Sxx)

# Use the middle column of the spectrogram to intialize the minimum and maximum values for the color map
middle_index =  len(times) // 2
middle_column = spec[:, middle_index]
vmin = 0  
vmax = np.max(middle_column) 

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

c = 311.1 # Default speed of sound, average of dataset, m/s
fa = np.max(coords_array[:, 1]) 
fr = np.min(coords_array[:, 1])
#insert method to get initial model here
fm = (fa+fr)/2 #- 20

#find the closest coordinate to f0
closest_index = np.argmin(np.abs(coords_array[:, 1] - fm))
f0 = coords_array[closest_index, 1] 
tprime0 = coords_array[closest_index, 0]  
t_hold = np.inf
for i,t in enumerate(coords_array[:, 0]):
    if t != tprime0:
        if (t - tprime0) < t_hold:
            t_hold = abs(t - tprime0)
            second_index = i

v0 = c*abs(fa-fr) / (2 * f0)
slope = (coords_array[closest_index,1] - coords_array[second_index,1]) / (coords_array[closest_index,0] - coords_array[second_index,0])
l = -((f0*v0**2/c)*(1-(v0/c)**2)**(-3/2))/slope 
m0 = [f0, v0, l, tprime0, c]

# Perform inversion using the initial model parameters and coordinates array
m0 = [f0, v0, l, tprime0, c]
sigma_prior = [40, 1, 1, 200, 1]
m,_,_, F_m = invert_f(m0,sigma_prior, coords_array, num_iterations=3)
m0[0] = m[0]
m0[3] = m[3]

sigma_f0 = 150
sigma_v0 = 100
sigma_l = 10000
sigma_tprime0 = 200
sigma_c = 100

m0 = [f0, v0, l, tprime0, c]
sigma_prior = [sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c]
m,_,_, F_m = invert_f(m0,[sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c], coords_array, num_iterations=3)
v0 = m[1]
l = m[2]
tprime0 = m[3]
c = m[4]

mprior = []
mprior.append(v0)
mprior.append(l)
mprior.append(tprime0)
mprior.append(c)

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


if len(peaks) <= 15:
    corridor_width = 10
else:
    corridor_width = 5

tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks,freqpeak, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax)
for o in range(len(f0_array)):
    mprior.append(float(f0_array[o]))

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
for p in range(len(f0_array)):
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

if abs(slope) < 1:
    sigma_prior = [10, 125, 15000, 30, 100]
else:
    sigma_prior = [10, 30, 500, 30, 100]

m, covm0, covm, f0_array, F_m = full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations=2, sigma=3, off_diagonal=False)
v0 = m[0]
l = m[1]
tprime0 = m[2]
c = m[3]
Cpost = np.sqrt(np.diag(covm))
Cpost0 = np.sqrt(np.diag(covm0))

# Find the index of the closest time value to tprime0
closest_index = np.argmin(np.abs(tprime0 - times))
# Get the arrival time values at the closest time index
arrive_time = spec[:, closest_index]
# Set negative arrival time values to 0
for i in range(len(arrive_time)):
    if arrive_time[i] < 0:
        arrive_time[i] = 0
# Plot settings and calculations
vmin = np.min(arrive_time) 
vmax = np.max(arrive_time)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     
ax1.plot(torg, data, 'k', linewidth=0.5)
ax1.set_title(title)

ax1.margins(x=0)
ax1.set_position([0.125, 0.6, 0.775, 0.3]) 
ax1.set_ylabel('Counts')
# Plot spectrogram
cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
ax2.set_xlabel('Time (s)')
f0lab = []
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
        for c_index  in range(4, len(Cpost0)):
            xmin = f0_array[c_index-4] - Cpost0[c_index]
            xmax = f0_array[c_index-4] + Cpost0[c_index]
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

if len(f0lab) > 10:
    # Split f0lab into lines of 10 entries each
    f0lab_lines = []
    for i in range(0, len(f0lab), 10):
        line = ', '.join(["%.2f" % f for f in f0lab[i:i+10]])
        f0lab_lines.append(line)
    f0lab_str = (',\n').join(f0lab_lines)
    f0lab_str = '[' + f0lab_str + ']'
else:
    f0lab_str = '[' + ', '.join(["%.2f" % f for f in f0lab]) + ']'

if isinstance(F_m, str):
        if med_df == "NaN":
            ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % Cpost0[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % Cpost0[0]+' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] + ' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u2080 = ' + f0lab_str + ' \u00B1 ' + "%.2f" % np.median(Cpost0[3:]) +' Hz,\n[' + F_m + ']', fontsize=fss)
        else:
            ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % Cpost0[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % Cpost0[0] +' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] + ' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u2080 = ' + f0lab_str + ' \u00B1 ' + "%.2f" % np.median(Cpost0[3:]) +' Hz, df\u2080 = ' + "%.2f" % med_df + ' \u00B1 ' + "%.2f" % mad_df + ' Hz\n[' + F_m + ']', fontsize=fss)
elif med_df == "NaN":
    ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % Cpost0[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % Cpost0[0]+' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] + ' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u2080 = ' + f0lab_str + ' \u00B1 ' + "%.2f" % np.median(Cpost0[3:]) +' Hz,\nMisfit: ' + "%.4f" % F_m, fontsize=fss)
else:
    ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % Cpost0[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % Cpost0[0] +' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] + ' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u2080 = ' + f0lab_str + ' \u00B1 ' + "%.2f" % np.median(Cpost0[3:]) +' Hz, df\u2080 = ' + "%.2f" % med_df + ' \u00B1 ' + "%.2f" % mad_df + ' Hz\nMisfit: ' + "%.4f" % F_m, fontsize=fss)

ax2.legend(loc='upper right',fontsize = 'small')
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
plt.show()
plt.close()

