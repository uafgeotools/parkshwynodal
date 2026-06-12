"""
This script calculates and visualizes results from changing variables in 
the Doppler effect equation. It shows the frequency signals perceived by a 
receiver from a moving source (aircraft). The varying parameters include 
velocity, altitude, wave propagation speed, time of closest approach, and 
emitted frequency. It generates subplots to show the impact of each parameter 
on the observed frequency over time.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Default doppler shift parameters
DEFAULTS = {
    'v': 80,
    'd0': 4000,
    't0': 120,
    'c': 320,
    'fs': 150
}

# Variables ranges for plotting differing results for changes in each variable
VAR_RANGES = {
    'base': np.arange(0, 241, 1),
    'v': np.arange(0, 200, 20),
    'd0': np.arange(0, 5000, 200),
    't0': np.arange(0, 240, 20),
    'c': np.arange(200, 500, 5),
    'fs': np.arange(0, 250, 20),
}

# Time array for plotting
time_receiver = VAR_RANGES['base']


def get_t(v, d0, t0, c, time_receiver, org=False):
    """Calculate time in aircraft reference frame from 
    spectrogram reference frame.

    :type v: numpy.int64 float
    :param v: Velocity of the aircraft
    :type d0: numpy.int64 float
    :param d0: Altitude of the aircraft
    :type c: numpy.int64 float
    :param c: Velocity of the wave propagation
    :type t0: numpy.int64 float
    :param t0: Time of closest approach in spectrogram reference frame
    :type time_receiver: np.array float
    :param time_receiver: Time array in spectrogram reference frame
    :type org: bool
    :param org: If True, return original time array
    :rtype: np.array
    :return: Time array in aircraft reference frame
    """
    beta = v/c
    if org:
        arg = time_receiver - t0
    else:
        arg = time_receiver - t0 + d0/c
    discriminant = arg**2 - (1 - beta**2) * (arg**2 - (d0/c)**2)
    return (arg - np.sqrt(discriminant)) / (1 - beta**2)


def get_f(v, d0, tprime_array, c, fs):
    """Calculate observed frequency (Doppler effect).

    :type v: numpy.int64 float
    :param v: Velocity of the aircraft
    :type d0: numpy.int64 float
    :param d0: Altitude of the aircraft
    :type c: numpy.int64 float
    :param c: Velocity of the wave propagation
    :type tprime_array: np.array
    :param tprime_array: Time array in aircraft reference frame
    :type fs: numpy.int64 float
    :param fs: Emitted frequency from aircraft
    :rtype: np.array
    :return: Frequency array of received frequency at seonsors corresponding 
    to tprime_array
    """
    base = np.sqrt(d0**2 + (v * tprime_array)**2)
    base_checked = np.where(base != 0, base, 1e-10)
    return fs / (1 + (v/c) * (v*tprime_array) / base_checked)


# Create subplots
fig, axs = plt.subplots(2, 3, figsize=(11, 7), sharex=True, sharey=True)
axs = axs.flatten()
cm = plt.cm.rainbow

# Plot for each variable
var_names = list(VAR_RANGES.keys())
for n, var_name in enumerate(var_names):
    params = DEFAULTS.copy()

    # Base plot with default parameters
    tprime_array = get_t(params['v'], params['d0'], params['t0'], params['c'], 
                   time_receiver)
    print(tprime_array)
    ft = get_f(params['v'], params['d0'], tprime_array, params['c'], 
               params['fs'])

    axs[n].plot(time_receiver, ft, c='k', linewidth=0.5, zorder=10)
    axs[n].axvline(params['t0'], c='k', linewidth=0.5, zorder=10)
    axs[n].set_title(f'Varying {var_name}')
    axs[n].set_ylim(100, 225)
    axs[n].set_xlim(0, 240)

    # Set labels
    if n in [0, 3]:
        axs[n].set_ylabel('Frequency (Hz)')
    if n in [3, 4, 5]:
        axs[n].set_xlabel('Time (s)')

    # Skip variable plotting for base case
    if n != 0:
        # Add colorbar for each subplot except the first
        var_values = VAR_RANGES[var_name]
        norm = plt.Normalize(var_values.min(), var_values.max())
        sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
        plt.colorbar(sm, ax=axs[n], orientation='vertical', pad=0.01, aspect=30)

        # Plot for varying parameters 
        for val in var_values:
            params[var_name] = val
            tprime_array = get_t(params['v'], params['d0'], params['t0'], 
                                 params['c'], time_receiver)
            ft = get_f(params['v'], params['d0'], tprime_array, params['c'], 
                    params['fs'])
            axs[n].plot(time_receiver, ft, color=cm(norm(val)), linewidth=0.5)
plt.tight_layout()
plt.show()