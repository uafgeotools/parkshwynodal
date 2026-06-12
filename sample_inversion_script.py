'''User can put any time series data from IRIS that contains a moving source  
and follow the prompts to perform a Doppler inversion. The script will guide the 
user through picking points on the spectrogram, associating those picks with 
overtone curves, and selecting a time window for the inversion. The final output 
will be a plot of the spectrogram with the inversion results overlaid, along 
with a title that includes the estimated model parameters and their 
uncertainties. The inversion results will be printed in the terminal for each 
iteration.'''

import sys
import numpy as np

from pathlib import Path
from obspy.core import UTCDateTime
from scipy.signal import spectrogram

# --- Fix sys.path ---
repo_root = Path(__file__).resolve().parents[0]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.get_save_data import load_waveform
from src.doppler_funcs import DopplerInversion as DI
from src.fig_func import SpecPlot, GetPicks

arrive_time = (UTCDateTime("2019-02-21T20:33:34") + 120).timestamp

data, sample_rate, t_wf, title = load_waveform('1011', arrive_time, 
                                               spec_window=120, component="Z")

# Compute spectrogram
WIN_LEN = 1  # window length, in s
NPER = int(WIN_LEN * sample_rate)
frequencies, times, Sxx = spectrogram(
    data, sample_rate, scaling='density', nperseg=NPER,
    noverlap=int(NPER * .9), detrend='constant'
)

spec_plot = SpecPlot(Sxx, sample_rate, t_wf, data, times, frequencies, 
                 arrive_time)
spec, MDF = spec_plot.remove_median()
middle_index = len(times) // 2
middle_column = spec[:, middle_index]
vmin, vmax = 0, np.max(middle_column)
base_dir = f'{repo_root}/input/data_picks/sample_script_picks/'

fpks1 = f'doppler_picks.txt'

spec_load = GetPicks(spec_plot, vmin, vmax, save_picks=True,review_picks=True)
coords_array, start_time = spec_load.single_doppler_data(fpks1, base_dir)

# Estimate initial model parameters from picked points
c = 320 # Speed of sound (m/s)
fa, fr = np.max(coords_array[:, 1]), np.min(coords_array[:, 1])
fm = (fa + fr) / 2
closest_index = np.argmin(np.abs(coords_array[:, 1] - fm))
fs, t0 = coords_array[closest_index, 1], coords_array[closest_index, 0]
t_hold, second_index = np.inf, None
for i, t in enumerate(coords_array[:, 0]):
    if t != t0 and abs(t - t0) < t_hold:
        t_hold = abs(t - t0)
        second_index = i
v = c * abs(fa - fr) / (2 * fs)  # Initial velocity estimate
slope = (
    (coords_array[closest_index, 1] - coords_array[second_index, 1]) /
    (coords_array[closest_index, 0] - coords_array[second_index, 0])
)
d0 = -(
    (fs * v ** 2 / c) * (1 - (v / c) ** 2) ** (-3 / 2)
) / slope  # Initial length estimate
m0 = [v, d0, t0, c, fs]
sigma_prior = [1, 1, 200, 1, 40]  # Initial prior uncertainties
fobs = []
tobs = []
for t, f in coords_array:
    tobs.append(t)
    fobs.append(f)

aircraft_inversion = DI(
    fobs, tobs, m0, sigma_prior, num_iterations=3, off_diagonal=False)
# First inversion to refine model
m, _, _, _, F_m = aircraft_inversion.full_inversion(
    [len(fobs)])
m0[4], m0[2] = m[4], m[2]

# Second inversion with wider priors
sigma_prior = [100, 10000, 200, 100, 150]

aircraft_inversion = DI(
    fobs, tobs, m0, sigma_prior, num_iterations=3, off_diagonal=False)

m, _, _, _, F_m = aircraft_inversion.full_inversion(
    [len(fobs)])
v, d0, t0, c = m[0], m[1], m[2], m[3]
mprior = [v, d0, t0, c]
fpks2 = f'overtone_picks.txt'
peaks, time_peak = spec_load.overtone_data(t0, file_name=fpks2, 
                                             BASE_DIR=base_dir)

# Automatically associate picked peaks with overtone curves
corridor_width = 10 if len(peaks) <= 15 else 5
tobs, fobs, peaks_assos, fs_array = spec_load.auto_picks_full(
    peaks, time_peak, corridor_width, t0, v, d0, c, sigma_prior
    )

mprior += [float(f) for f in fs_array]
fpks3 = f'time_picks.txt'
tobs, fobs, peaks_assos = spec_load.final_data(tobs, fobs, fs_array, 
                                               peaks_assos, file_name=fpks3, 
                                               BASE_DIR=base_dir)

# Final inversion using filtered picks
sigma_prior = (
    [10, 125, 15000, 30, 100] if abs(slope) < 1
    else [10, 30, 500, 30, 100]
)

aircraft_inversion = DI(fobs, tobs, mprior, sigma_prior, 
                        num_iterations=4, off_diagonal=False)

m, covm0, covm, fs_array, F_m = aircraft_inversion.full_inversion(
    peaks_assos, sigma=3)

v, d0, t0, c = m[0], m[1], m[2], m[3]
Cpost, Cpost0 = np.sqrt(np.diag(covm)), np.sqrt(np.diag(covm0))

spec_load.plot_spectrogram(
    title, t0, v, d0, c, fs_array, F_m, Cpost0, middle_index, file_name=None, 
    plot_show=True, gt = True)