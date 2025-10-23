import obspy
import sys
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timezone
from scipy.signal import spectrogram
from pathlib import Path

# Ensure repository root is on sys.path so local package 'src' can be imported
repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
from src.doppler_funcs import calc_ft, invert_f
from src.main_inv_fig_functions import  remove_median

wave_form_data_path = "/scratch/naalexeev/NODAL/"
c = 320 #Speed of sound in m/s
start_time = 1550158642.26246 #Start time of spectrogram

#Extract time data 
ht = datetime.fromtimestamp(start_time, tz=timezone.utc)                      
h = ht.hour
mins = ht.minute
secs = ht.second
month = ht.month
day = ht.day
h_u = str(h+1)

#load waveform data
p = wave_form_data_path + "/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day)+"T"+str(h_u)+":00:00.000000Z.1173.mseed"
tr = obspy.read(p)
tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs , tr[2].stats.starttime + (mins * 60) + secs + 240)
data = tr[2][:]
fs = int(tr[2].stats.sampling_rate)
title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
t_wf = tr[2].times()

# Compute spectrogram
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 

# Remove median background spectrum calulated for entire 340 sec window
spec, MDF = remove_median(Sxx)

#Get vmin and vmax of the center of the signal for consistent color scale
middle_index =  len(times) // 2
middle_column = spec[:, middle_index]
vmin = 0  
vmax = np.max(middle_column) 

#Prepicked t' and F values to get prior model
x = [112.48911983478979, 59.65932080234049, 186.52395930932946, 102.98341205040444, 120.34960896418536]
y = [140.02964002964, 188.29218829218826, 93.7170937170937, 153.9234039234039, 128.81712881712878]

coords = [(x[i], y[i]) for i in range(len(x))]
coords_array = np.array(coords)

# Create a subplot for the visualization
fig, ax = plt.subplots(5,1,figsize=(8/1.4, 14/1.4),sharex=False)
cax = ax[0].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[0].axhline(y=coords_array[1,1], color='black', linestyle='--', linewidth=1)
ax[0].axhline(y=coords_array[2,1], color='black', linestyle='--', linewidth=1)
ax[0].axhline(y=(coords_array[1,1]+coords_array[2,1])/2, color='red', linestyle='--', linewidth=0.7)
ax[0].axvline(x=coords_array[0,0], color='red', linestyle='--', linewidth=0.7)
slope = (coords_array[4,1] - coords_array[3,1]) / (coords_array[4,0] - coords_array[3,0])

# Add points at x=70 and x=150 using the slope for line extension
y_70 = coords_array[3,1] + slope * (70 - coords_array[3,0])
y_150 = coords_array[3,1] + slope * (150 - coords_array[3,0])
ax[0].plot([70, 150], [y_70, y_150], color='blue', linestyle='--', linewidth=1, zorder=1)

# Move scatter plots after all lines so they appear on top
ax[0].scatter(coords_array[1:3, 0], coords_array[1:3, 1], c='black', marker='x', s=100, linewidths=3,label="f_initial + f_final")
ax[0].scatter(coords_array[0, 0], coords_array[0, 1], c='red', marker='x', s=100, linewidths=3, label="t'0 + f0")
ax[0].scatter(coords_array[3:5, 0], coords_array[3:5, 1], c='blue', marker='x', s=100, linewidths=3, label="Slope of l")
ax[0].set_ylabel('Frequency (Hz)')
ax[0].set_title("(a) data picks to get prior model", fontsize='small')

cax = ax[1].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

#Method to get initial model 
f0 = ((coords_array[1,1]+coords_array[2,1])/2) * 0.84
t0 = coords_array[0,0] 
del_f = coords_array[1,1]-coords_array[2,1]
v0 = (c/del_f)*(np.sqrt((f0**2+del_f**2)) - f0) 
slope_t0prime = slope*((1-(v0/c)**2)**(-3/2))
l = -(f0*v0**2/(c*slope_t0prime))

#Plot initial model
m0 = [f0, v0, l, t0, c]
print('Initial model:', m0)
ft = calc_ft(times, m0[3], m0[0], m0[1], m0[2], m0[4])
ax[1].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[1].scatter(coords_array[1:3, 0], coords_array[1:3, 1], c='black', marker='x', s=100, linewidths=3,label="f_initial + f_final")
ax[1].scatter(coords_array[0, 0], coords_array[0, 1], c='red', marker='x', s=100, linewidths=3, label="t'0 + f0")
ax[1].scatter(coords_array[3:5, 0], coords_array[3:5, 1], c='blue', marker='x', s=100, linewidths=3, label="Slope of l")
ax[1].set_ylabel('Frequency (Hz)')
ax[1].set_title("(b) prior model", fontsize='small')

m0 = [f0, v0, l, t0, c]
sigma_prior = [40, 1, 1, 200, 1]
m,_,_, F_m = invert_f(m0,sigma_prior, coords_array, num_iterations=3)
m0[0] = m[0]
m0[3] = m[3]

tf = np.arange(0, 240, 1)

sigma_f0 = 150
sigma_v0 = 100
sigma_l = 10000
sigma_t0 = 200
sigma_c = 100
sigma_prior = [sigma_f0, sigma_v0, sigma_l, sigma_t0, sigma_c]
# First inversion using all 5 picked data
m, covm,_, F_m = invert_f(m0, sigma_prior, coords_array, num_iterations=5)
ft = calc_ft(times, m[3], m[0], m[1], m[2], m[4])
#Plot model derived from all 5 picked data
cax = ax[2].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[2].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[2].set_ylabel('Frequency (Hz)')
ax[2].set_title("(c) measured model => prior model", fontsize='small')

#Create corridor around prior model and extract data within corridor
peaks = []
coord_inv = []
upper_array = []
lower_array = []
corridor_width = 10 
time_corr = np.arange(0, 240, 1)
for ttt in time_corr:
    t_f = (np.abs(times - ttt)).argmin()
    upper = int(ft[t_f] + corridor_width)
    lower = int(ft[t_f] - corridor_width)
    if lower < 0:
        lower = 0
    if upper > len(frequencies):
        upper = len(frequencies)
    tt = spec[lower:upper, t_f]
    try:
        max_amplitude_index = np.argmax(tt)
    except:
        continue
    max_amplitude_frequency = frequencies[max_amplitude_index+lower]
    peaks.append(max_amplitude_frequency)
    coord_inv.append((times[t_f], max_amplitude_frequency))
    upper_array.append(upper)
    lower_array.append(lower)

coord_inv_array = np.array(coord_inv)
cax = ax[3].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[3].plot(coord_inv_array[:, 0], np.array(upper_array), 'r', linewidth=1)
ax[3].plot(coord_inv_array[:, 0], np.array(lower_array), 'r', linewidth=1)
ax[3].set_title("(d) data extracted from model corridor (prior model \u00B1 10)", fontsize='small')

#With all data extracted from corridor, do inversion and remove outliers
prior_sigma = [10,30,600,30,80] #prior sigma values for f0, v0, l, t0, c
m,_,_,F_m = invert_f(m, prior_sigma, coord_inv_array, num_iterations=3)

f0 = m[0]
v0 = m[1]
l = m[2]
t0 = m[3]
c = m[4]
ft = calc_ft(time_corr, m[3], m[0], m[1], m[2], m[4])

delf = np.array(ft) - np.array(peaks)

new_coord_inv_array = []
for i in range(len(delf)):
    if np.abs(delf[i]) <= 3:
        new_coord_inv_array.append(coord_inv_array[i])
coord_inv_array = np.array(new_coord_inv_array)

ax[3].scatter(coord_inv_array[:, 0], coord_inv_array[:, 1], c='black', marker='x', s=20)
ax[3].set_ylabel('Frequency (Hz)')

prior_sigma = [10,30,500,30,80] #prior sigma values for f0, v0, l, t0, c
m,covm0,covm_norm,F_m = invert_f(m, prior_sigma, coord_inv_array, num_iterations=6, sigma=2)
covm0 = np.sqrt(np.diag(covm0))

f0 = m[0]
v0 = m[1]
l = m[2]
t0 = m[3]
c = m[4]

ft = calc_ft(times, t0, f0, v0, l, c)
cax = ax[4].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[4].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[4].set_ylabel('Frequency (Hz)')
ax[4].set_title("(e) posterior model", fontsize='small')

#make all axis tick labels smaller
for i in range(5):
	ax[i].tick_params(axis='both', which='major', labelsize='x-small')
	ax[i].tick_params(axis='both', which='minor', labelsize='x-small')

#make  the gap between subplots smaller
plt.subplots_adjust(hspace=0.3)
ax[4].set_xlabel('Time (s)')
plt.tight_layout()

fig.savefig("inversion_steps.jpg", dpi=600)
plt.close()
