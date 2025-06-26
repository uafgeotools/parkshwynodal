import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj
from prelude import *
from scipy.optimize import minimize

# Function to calculate L1 norm and produce best fit line for plot
def fit_l1_line(x, y, bounds=None):
    """
    Fits a line to the given data points (x, y) using the L1 norm (minimizing absolute deviations).

    Parameters:
        x (list or np.array): Independent variable data points.
        y (list or np.array): Dependent variable data points.
        bounds (tuple, optional): A tuple (min_val, max_val) specifying the range of values to consider for both x and y.

    Returns:
        tuple: Slope (m) and intercept (b) of the best fit line.
    """
    # Apply bounds if provided
    if bounds is not None:
        min_val, max_val = bounds
        x = np.array(x)  # Ensure x is a NumPy array
        y = np.array(y)  # Ensure y is a NumPy array
        mask = (x >= min_val) & (x <= max_val) & (y >= min_val) & (y <= max_val)
        x = np.array(x)[mask]
        y = np.array(y)[mask]

    # Define the objective function for L1 norm
    def L1(params):
        m, b = params
        x_array = np.array(x)  # Ensure x is a NumPy array
        return np.sum(np.abs(y - (m * x_array + b)))

    # Initial guess for slope and intercept
    initial_guess = [0, 0]

    # Minimize the L1 norm
    result = minimize(L1, initial_guess)

    # Extract the slope and intercept from the result
    m, b = result.x
    return m, b

file = open('output/DH8Adata_atmosphere_full.txt','r')

inverse_times = []
inverse_dists = []
inverse_speeds = []
flight_nums = []
comp_times = []
c_array = []
error_vel = []
error_dist = []
error_time = []

for line in file.readlines():
    lines = line.split(',')
    inverse_times.append(float(lines[4]))
    
    inverse_dists.append(abs(float(lines[6])))
    inverse_speeds.append(float(lines[5]))
    comp_times.append(float(lines[3])) 
    flight_num = lines[1]
    c_array.append(float(lines[11]))
    # Process old and new peaks
    error = np.array(lines[8])
    error = str(error)
    error = np.char.replace(error, '[', '')
    error = np.char.replace(error, ']', '')
    error = str(error)
    error = np.array(error.split(' '))

    error_vel.append(float(error[0]))
    error_dist.append(float(error[1]))
    error_time.append(float(error[2]))


file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

fr_times = []
fr_dists = []
fr_speeds = []

for line in file_in.readlines():
    text = line.split(',')
    flight_id = text[1]
    closest_time = float(text[5])
    if flight_id not in flight_nums and closest_time not in comp_times:
        continue
    index = comp_times.index(closest_time)
    c = c_array[index]
    fr_dists.append(abs(np.sqrt(float(text[4])**2 + float(text[6])**2)))

    fr_speeds.append(float(text[7]))
    ta_old = calc_time(float(closest_time),float(text[4]),float(text[5]),float(c))
    fr_times.append(120)

fig, axs = plt.subplots(1, 3, figsize=(24, 6), sharey=False, layout='constrained')
scatter1 = axs[0].scatter(inverse_speeds, fr_speeds, c='k', s=15, zorder=2)
axs[0].errorbar(inverse_speeds, fr_speeds, xerr=error_vel, fmt='none', c='k', zorder=1)
axs[0].set_title("Velocity (m/s)", fontsize=10)
axs[0].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0].set_aspect('equal')
axs[0].set_xlabel('Nodal Data', fontsize=8)
axs[0].set_ylabel('Flightradar24', fontsize=8)
axs[0].tick_params(axis='both', labelsize=8)
m, b = fit_l1_line(inverse_speeds, fr_speeds)
x = np.linspace(min(inverse_speeds), max(inverse_speeds), 100)
axs[0].plot(x, m * x + b, color='k')

scatter2 = axs[1].scatter(inverse_dists, fr_dists, c='k', s=15, zorder=2)
axs[1].errorbar(inverse_dists, fr_dists, xerr=error_dist, fmt='none', c='k', zorder=1)
axs[1].set_title("Distance (m)", fontsize=10)
axs[1].axline((0, 0), slope=1, color='black', linestyle='--')
axs[1].set_aspect('equal', adjustable='box')
axs[1].set_xlabel('Nodal Data', fontsize=8)
axs[1].set_ylabel('Flightradar24', fontsize=8)
axs[1].tick_params(axis='both', labelsize=8)
m, b = fit_l1_line(inverse_dists, fr_dists)
x = np.linspace(min(inverse_dists), max(inverse_dists), 100)
axs[1].plot(x, m * x + b, color='k')

scatter3 = axs[2].scatter(inverse_times, fr_times, c='k', s=15, zorder=2)
axs[2].errorbar(inverse_times, fr_times, xerr=error_time, fmt='none', c='k', zorder=1)
axs[2].set_title("Time (s)", fontsize=10)
axs[2].set_xlabel('Nodal Data', fontsize=8)
axs[2].set_ylabel('Flightradar24', fontsize=8)
axs[2].tick_params(axis='both', labelsize=8)

for i, ax in enumerate(axs):
    ax.set_box_aspect(1)  # Make Time plots square
plt.show()
