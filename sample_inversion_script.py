import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from src.main_inv_fig_functions import remove_median, get_auto_picks_full
from src.doppler_funcs import invert_f, calc_ft, full_inversion
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from matplotlib.ticker import MaxNLocator

# Interactive picking of points on spectrogram for overtone curve
def pick_points_on_spectrogram(times, frequencies, spec, vmin, vmax, prompt, axvline=None):
    coords = []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    if axvline is not None:
        plt.axvline(x=axvline, c='#377eb8', ls='--')
    def onclick(event):
        if event.xdata is not None and event.ydata is not None:
            coords.append((event.xdata, event.ydata))
            plt.scatter(event.xdata, event.ydata, color='black', marker='x')
            plt.draw()
            print('Clicked:', event.xdata, event.ydata)
    plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show(block=True)
    return coords

# Interactive picking of single points (overtone peaks)
def pick_single_points(times, frequencies, spec, vmin, vmax, prompt, axvline=None):
    peaks, freqpeak = [], []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    if axvline is not None:
        plt.axvline(x=axvline, c='#377eb8', ls='--')
    def onclick(event):
        if event.xdata is not None and event.ydata is not None:
            peaks.append(event.ydata)
            freqpeak.append(event.xdata)
            plt.scatter(event.xdata, event.ydata, color='black', marker='x')
            plt.draw()
            print('Clicked:', event.xdata, event.ydata)
    plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show(block=True)
    return peaks, freqpeak

# Interactive picking of time window for inversion
def pick_time_window(times, frequencies, spec, vmin, vmax, tobs, fobs):
    set_time = []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    plt.scatter(tobs, fobs, color='black', marker='x')
    def onclick(event):
        if event.xdata is not None:
            set_time.append(event.xdata)
            plt.scatter(event.xdata, event.ydata, color='red', marker='x')
            plt.draw()
            print('Clicked:', event.xdata, event.ydata)
    plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show(block=True)
    return set_time

# Download waveform data from IRIS PH5WS
client = Client("http://service.iris.edu", service_mappings={"dataselect": "http://service.iris.edu/ph5ws/dataselect/1"})
starttime = UTCDateTime("2019-02-21T20:33:34")
endtime = UTCDateTime("2019-02-21T20:37:34")
st = client.get_waveforms("ZE", "1011", "*", "DPZ", starttime, endtime)
tr = st[0]
data = tr.data
t_wf = tr.times()
fs = int(tr.stats.sampling_rate)
title = f'{tr.stats.network}.{tr.stats.station}.{tr.stats.location}.{tr.stats.channel} − starting {tr.stats["starttime"]}'

# Compute spectrogram
WIN_LEN = 1  # window length, in s
NPER = int(WIN_LEN * fs)
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=NPER, noverlap=int(NPER * .9), detrend='constant')

spec, MDF = remove_median(Sxx)  # Remove median for better visualization
middle_index = len(times) // 2
middle_column = spec[:, middle_index]
vmin, vmax = 0, np.max(middle_column)

# User picks overtone curve points
print("Please pick the points on the spectrogram that correspond to the primary overtone of the doppler curves.")
while True:
    coords = pick_points_on_spectrogram(times, frequencies, spec, vmin, vmax, "Pick overtone curve points")
    if input("Do you want to repick your points? (y or n)").lower() != 'y':
        break
coords_array = np.array(coords)
print(coords_array)


# Estimate initial model parameters from picked points
c = 320 #11.1  # Speed of sound (m/s)
fa, fr = np.max(coords_array[:, 1]), np.min(coords_array[:, 1])  # Max/min frequency
fm = (fa + fr) / 2
closest_index = np.argmin(np.abs(coords_array[:, 1] - fm))
f0, t0 = coords_array[closest_index, 1], coords_array[closest_index, 0]
t_hold, second_index = np.inf, None
for i, t in enumerate(coords_array[:, 0]):
    if t != t0 and abs(t - t0) < t_hold:
        t_hold = abs(t - t0)
        second_index = i
v0 = c * abs(fa - fr) / (2 * f0)  # Initial velocity estimate
slope = (coords_array[closest_index, 1] - coords_array[second_index, 1]) / (coords_array[closest_index, 0] - coords_array[second_index, 0])
l = -((f0 * v0 ** 2 / c) * (1 - (v0 / c) ** 2) ** (-3 / 2)) / slope  # Initial length estimate
m0 = [f0, v0, l, t0, c]
sigma_prior = [40, 1, 1, 200, 1]  # Initial prior uncertainties

# First inversion to refine model
m, _, _, F_m = invert_f(m0, sigma_prior, coords_array, num_iterations=3)
m0[0], m0[3] = m[0], m[3]

# Second inversion with wider priors
sigma_prior = [150, 100, 10000, 200, 100]
m, _, _, F_m = invert_f(m0, sigma_prior, coords_array, num_iterations=3)
v0, l, t0, c = m[1], m[2], m[3], m[4]
mprior = [v0, l, t0, c]

# User picks overtone peaks
print("Please pick one point on each overtone, it does not have to be at the center of the doppler.")
while True:
    peaks, freqpeak = pick_single_points(times, frequencies, spec, vmin, vmax, "Pick overtone peaks", axvline=t0)
    if input("Do you want to repick your points? (y or n)").lower() != 'y':
        break

# Automatically associate picked peaks with overtone curves
corridor_width = 10 if len(peaks) <= 15 else 5
tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks, freqpeak, times, frequencies, spec, corridor_width, t0, v0, l, c, sigma_prior, vmax)
mprior += [float(f) for f in f0_array]

# User picks time window for inversion
print('Please pick two points on the spectrogram that correspond to the start and end of the time window you want pull data from in the inversion.')
while True:
    set_time = pick_time_window(times, frequencies, spec, vmin, vmax, tobs, fobs)
    if input("Do you want to repick your points? (y or n)").lower() != 'y':
        break
start_time, end_time = set_time[:2]

# Filter picks to only those within the selected time window
ftobs, ffobs, peak_ass = [], [], []
cum = 0
for p in range(len(f0_array)):
    count = 0
    for j in range(cum, cum + peaks_assos[p]):
        if start_time <= tobs[j] <= end_time:
            ftobs.append(tobs[j])
            ffobs.append(fobs[j])
            count += 1
    cum += peaks_assos[p]
    peak_ass.append(count)
peaks_assos = peak_ass
tobs, fobs = ftobs, ffobs

# Final inversion using filtered picks
sigma_prior = [10, 125, 15000, 30, 100] if abs(slope) < 1 else [10, 30, 500, 30, 100]
m, covm0, covm, f0_array, F_m = full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations=2, sigma=3, off_diagonal=False)
v0, l, t0, c = m[0], m[1], m[2], m[3]
Cpost, Cpost0 = np.sqrt(np.diag(covm)), np.sqrt(np.diag(covm0))

# Plot results
closest_index = np.argmin(np.abs(t0 - times))
arrive_time = np.clip(spec[:, closest_index], 0, None)
vmin, vmax = np.min(arrive_time), np.max(arrive_time)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8, 6))

# Plot raw waveform
ax1.plot(t_wf, data, 'k', linewidth=0.5)
ax1.set_title(title)
ax1.margins(x=0)
ax1.set_position([0.125, 0.6, 0.775, 0.3])
ax1.set_ylabel('Counts')

# Plot spectrogram and inversion results
cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax2.set_xlabel('Time (s)')
t0prime = t0 + l/c
ax2.axvline(x=t0prime, c = '#377eb8', ls = '--',linewidth=0.5,label= "t\u2080' = " + "%.2f" % t0prime +' s')
ax2.axvline(x=t0, c = '#e41a1c', ls = '--', linewidth=0.7,label= "t\u2080 = " + "%.2f" % t0 +' s')
f0lab = sorted(f0_array)
for pp in range(len(f0_array)):
    f0 = f0_array[pp]
    ft = calc_ft(times, t0, f0, v0, l, c)

    ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7)
    ax2.scatter(t0prime, f0, color='black', marker='x', s=30, zorder=10)
fss = 'x-small'

# Estimate overtone frequency spacing and uncertainty
if len(f0lab) > 1:
    f_range = []
    NTRY = 1000
    for _ in range(NTRY):
        ftry = [f0_array[i - 4] + np.random.uniform(-Cpost0[i], Cpost0[i]) for i in range(4, len(Cpost0))]
        ftry = np.sort(ftry)
        f1 = [ftry[g] - ftry[g - 1] for g in range(1, len(ftry))]
        f_range.append(np.nanmedian(f1))
    med_df = np.nanmedian(f_range)
    mad_df = np.nanmedian(np.abs(f_range - med_df))
else:
    med_df = mad_df = "NaN"

# Format overtone frequencies for display
if len(f0lab) > 10:
    f0lab_lines = [', '.join([f"{f:.2f}" for f in f0lab[i:i + 10]]) for i in range(0, len(f0lab), 10)]
    f0lab_str = '[%s]' % (',\n'.join(f0lab_lines))
else:
    f0lab_str = '[' + ', '.join([f"{f:.2f}" for f in f0lab]) + ']'

# Compose plot title with inversion results and uncertainties
if isinstance(F_m, str):
    misfit_str = f"\n[{F_m}]"
else:
    misfit_str = f"\nMisfit: {F_m:.4f}"
df_str = f", df\u2080 = {med_df:.2f} \u00B1 {mad_df:.2f} Hz" if med_df != "NaN" else ""
ax2.set_title(
    f"t\u2080'= {t0:.2f} \u00B1 {Cpost0[2]:.2f} s, v\u2080 = {v0:.2f} \u00B1 {Cpost0[0]:.2f} m/s, "
    f"c = {c:.2f} \u00B1 {Cpost0[3]:.2f} m/s, l = {l:.2f} \u00B1 {Cpost0[1]:.2f} m, \n"
    f"f\u2080 = {f0lab_str} \u00B1 {np.median(Cpost0[3:]):.2f} Hz{df_str}{misfit_str}",
    fontsize=fss
)
ax2.legend(loc='upper right', fontsize='small')
ax2.set_ylabel('Frequency (Hz)')
ax2.margins(x=0)

# Add colorbar for spectrogram
ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])
cbar = plt.colorbar(mappable=cax, cax=ax3)
cbar.locator = MaxNLocator(integer=True)
cbar.update_ticks()
ax3.set_ylabel('Relative Amplitude (dB)')
ax2.margins(x=0)
#ax2.set_xlim(0, 240)
ax2.set_ylim(0, int(fs / 2))
ax1.tick_params(axis='both', which='major', labelsize=9)
ax2.tick_params(axis='both', which='major', labelsize=9)
ax3.tick_params(axis='both', which='major', labelsize=9)

# Overlay median-detrended frequency (MDF) for reference
spec2 = 10 * np.log10(MDF)
middle_column2 = spec2[:, middle_index]
vmin2, vmax2 = np.min(middle_column2), np.max(middle_column2)
ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2)
ax4.plot(middle_column2, frequencies, c='#ff7f00')
ax4.set_ylim(0, int(fs / 2))
ax4.set_xlim(vmax2 * 1.1, vmin2)
ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
ax4.grid(axis='y')
plt.show()
plt.close()
