import numpy as np
import sys
from pathlib import Path
from matplotlib import pyplot as plt
from scipy.signal import spectrogram


# --- Fix sys.path ---
repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.get_save_data import load_waveform
from src.doppler_funcs import DopplerInversion, DopplerCalc
from src.fig_func import SpecPlot

STATION = 1173
c = 320
crossing_time = 1550158642.26246 + 120

data, sample_rate, t_wf, title = load_waveform(
    STATION, crossing_time, spec_window=120
    )

# Compute spectrogram
frequencies, times, Sxx = spectrogram(
    data, sample_rate, scaling='density', nperseg=sample_rate, 
    noverlap=sample_rate * .9, detrend = 'constant'
    ) 

spec_plot = SpecPlot(
    Sxx, sample_rate, t_wf, data, times, frequencies, crossing_time
    )

spec, _ = spec_plot.remove_median()

middle_index =  len(times) // 2
middle_column = spec[:, middle_index]
vmin = 0  
vmax = np.max(middle_column) 

x = [112.48911983478979, 59.65932080234049, 186.52395930932946, 
     102.98341205040444, 120.34960896418536]
y = [140.02964002964, 188.29218829218826, 93.7170937170937, 
     153.9234039234039, 128.81712881712878]

coords = [(x[i], y[i]) for i in range(len(x))]
coords_array = np.array(coords)
fig_num = 5

# Create a subplot for the visualization
fig, ax = plt.subplots(fig_num,1,figsize=(8/1.4, 14/1.4),sharex=False)
cax = ax[0].pcolormesh(
    times, frequencies, spec, shading='gouraud', cmap='pink_r', 
    vmin=vmin, vmax=vmax
    )
ax[0].axhline(y=coords_array[1,1], color='black', linestyle='--', linewidth=1)
ax[0].axhline(y=coords_array[2,1], color='black', linestyle='--', linewidth=1)
ax[0].axhline(
    y=(coords_array[1,1]+coords_array[2,1])/2, color='red', linestyle='--', 
    linewidth=0.7
    )
ax[0].axvline(x=coords_array[0,0], color='red', linestyle='--', linewidth=0.7)
slope = (coords_array[4,1] 
         - coords_array[3,1]) / (coords_array[4,0] 
                                 - coords_array[3,0])

# Add points at x=70 and x=150 using the slope
y_70 = coords_array[3,1] + slope * (70 - coords_array[3,0])
y_150 = coords_array[3,1] + slope * (150 - coords_array[3,0])
ax[0].plot(
    [70, 150], [y_70, y_150], color='blue', linestyle='--', linewidth=1, 
    zorder=1
    )

# Move scatter plots after all lines so they appear on top
ax[0].scatter(
    coords_array[1:3, 0], coords_array[1:3, 1], c='black', marker='x', s=100, 
    linewidths=3,label="f_initial + f_final"
    )
ax[0].scatter(
    coords_array[0, 0], coords_array[0, 1], c='red', marker='x', s=100, 
    linewidths=3, label="t'0 + f0"
    )
ax[0].scatter(
    coords_array[3:5, 0], coords_array[3:5, 1], c='blue', marker='x', s=100, 
    linewidths=3, label="Slope of l"
    )
ax[0].set_ylabel('Frequency (Hz)')
ax[0].set_title("(a) data picks to get prior model", fontsize='small')

cax = ax[1].pcolormesh(
    times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, 
    vmax=vmax
    )
#insert method to get initial model here
fs = ((coords_array[1,1]+coords_array[2,1])/2) * 0.84
t0 = coords_array[0,0] 
del_f = coords_array[1,1]-coords_array[2,1]
v = (c/del_f)*(np.sqrt((fs**2+del_f**2)) - fs) 
slope_t0prime = slope*((1-(v/c)**2)**(-3/2))
d0 = -(fs*v**2/(c*slope_t0prime))

m0 = [v, d0, t0, c, fs]
print('Initial model:', m0)
doppler_values = DopplerCalc(v, d0, t0, c)
ft = doppler_values.calc_ft(times, fs)
ax[1].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[1].scatter(
    coords_array[1:3, 0], coords_array[1:3, 1], c='black', marker='x', s=100, 
    linewidths=3,label="f_initial + f_final"
    )
ax[1].scatter(
    coords_array[0, 0], coords_array[0, 1], c='red', marker='x', s=100, 
    linewidths=3, label="t'0 + f0"
    )
ax[1].scatter(
    coords_array[3:5, 0], coords_array[3:5, 1], c='blue', marker='x', s=100, 
    linewidths=3, label="Slope of l"
    )
ax[1].set_ylabel('Frequency (Hz)')
ax[1].set_title("(b) prior model", fontsize='small')

m0 = [v, d0, t0, c, fs]
sigma_prior = [40, 1, 1, 200, 1]
fobs = []
tobs = []
for t, f in coords_array:
    tobs.append(t)
    fobs.append(f)
aircraft_inversion = DopplerInversion(
    fobs, tobs, m0, sigma_prior, num_iterations=3, off_diagonal=False)
    
# First inversion to refine model
m, _, _, _, F_m = aircraft_inversion.full_inversion([len(fobs)])

m0[4] = m[4]
m0[2] = m[2]


tf = np.arange(0, 240, 1)

sigma_fs = 150
sigma_v = 100
sigma_d0 = 10000
sigma_t0 = 200
sigma_c = 100
sigma_prior = [sigma_v, sigma_d0, sigma_t0, sigma_c, sigma_fs]
aircraft_inversion = DopplerInversion(
    fobs, tobs, m0, sigma_prior, num_iterations=3, off_diagonal=False)
    
# First inversion to refine model
m, _, _, _, F_m = aircraft_inversion.full_inversion([len(fobs)])

doppler_values = DopplerCalc(m[0], m[1], m[2], m[3])
ft = doppler_values.calc_ft(times, m[4])

cax = ax[2].pcolormesh(
    times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, 
    vmax=vmax
    )
ax[2].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[2].set_ylabel('Frequency (Hz)')
ax[2].set_title("(c) measured model => prior model", fontsize='small')

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
cax = ax[3].pcolormesh(
    times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, 
    vmax=vmax
    )

ax[3].plot(coord_inv_array[:, 0], np.array(upper_array), 'r', linewidth=1)
ax[3].plot(coord_inv_array[:, 0], np.array(lower_array), 'r', linewidth=1)
ax[3].set_title(
    "(d) data extracted from model corridor (prior model \u00B1 10)", 
    fontsize='small'
    )

fobs = []
tobs = []
for t, f in coord_inv_array:
    tobs.append(t)
    fobs.append(f)
sigma_prior = [30, 600, 30, 80, 10] #prior sigma values for v, d0, t0, c, fs
aircraft_inversion = DopplerInversion(
    fobs, tobs, m, sigma_prior, num_iterations=3, off_diagonal=False)
    
# First inversion to refine model
m, _, _, _, F_m = aircraft_inversion.full_inversion([len(fobs)])

v = m[0]
d0 = m[1]
t0 = m[2]
c = m[3]
fs = m[4]

doppler_values = DopplerCalc(v, d0, t0, c)
ft = doppler_values.calc_ft(time_corr, fs)

delf = np.array(ft) - np.array(peaks)

new_coord_inv_array = []
for i in range(len(delf)):
    if np.abs(delf[i]) <= 3:
        new_coord_inv_array.append(coord_inv_array[i])
coord_inv_array = np.array(new_coord_inv_array)

ax[3].scatter(
    coord_inv_array[:, 0], coord_inv_array[:, 1], c='black', marker='x', s=20
    )
ax[3].set_ylabel('Frequency (Hz)')

fobs = []
tobs = []
for t, f in coord_inv_array:
    tobs.append(t)
    fobs.append(f)
sigma_prior = [30, 500, 30, 80, 10] #prior sigma values for v, d0, t0, c, fs
aircraft_inversion = DopplerInversion(
    fobs, tobs, m, sigma_prior, num_iterations=3, off_diagonal=False)
    
# First inversion to refine model
m, covm0, covm, fs_array, F_m = aircraft_inversion.full_inversion([len(fobs)])

covm0 = np.sqrt(np.diag(covm0))

v = m[0]
d0 = m[1]
t0 = m[2]
c = m[3]
fs = m[4]

doppler_values = DopplerCalc(v, d0, t0, c)
ft = doppler_values.calc_ft(times, fs)
cax = ax[4].pcolormesh(
    times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin,
    vmax=vmax
    )
ax[4].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[4].set_ylabel('Frequency (Hz)')
ax[4].set_title("(e) posterior model", fontsize='small')

#make all axis tick labels smaller
for i in range(fig_num):
	ax[i].tick_params(axis='both', which='major', labelsize='x-small')
	ax[i].tick_params(axis='both', which='minor', labelsize='x-small')
#make  the gap between subplots smaller
plt.subplots_adjust(hspace=0.3)
ax[fig_num-1].set_xlabel('Time (s)')
plt.tight_layout()
plt.show()
fig.savefig("inversion_steps.jpg", dpi=600)
plt.close()
