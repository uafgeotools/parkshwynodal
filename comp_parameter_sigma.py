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
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
stations = seismo_data['Station']
elevations = seismo_data['Elevation']
file_names = ['50_100_1000_30_1.txt','50_50_600_30_80.txt','50_10_200_30_80.txt','50_1_1_30_160.txt']

fr_times = []
fr_dists = []
fr_speeds = []
cc_array = []
error_bar = False
absolute_time = False
fig, axs = plt.subplots(4, 4, figsize=(15, 15), sharey=False, layout='constrained')
for ii, ff in enumerate(file_names):
    file = open(ff,'r')
    inverse_times = []
    inverse_dists = []
    inverse_speeds = []
    flight_nums = []
    comp_times = []
    c_array = []
    if error_bar:
        error_vel = []
        error_dist = []
        error_time = []
        error_c = []
    for line in file.readlines():
        lines = line.split(',')
        if absolute_time:
            inverse_times.append(float(lines[7]))
        else:
            inverse_times.append(float(lines[6]))
        inverse_dists.append(abs(float(lines[5])))
        inverse_speeds.append(abs(float(lines[4])))
        comp_times.append(float(lines[3])) 
        flight_num = lines[1]
        c_array.append(float(lines[8]))
        if error_bar:
            error = np.array(lines[-6])
            error = str(error)
            error = np.char.replace(error, '[', '')
            error = np.char.replace(error, ']', '')
            error = str(error)
            error = np.array(error.split(' '))

            error_vel.append(float(error[0]))
            error_dist.append(float(error[1]))
            error_time.append(float(error[2]))
            error_c.append(float(error[3]))

    if ii == 0:
        file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

        for line in file_in.readlines():
            text = line.split(',')
            date = text[0]

            month = int(date[4:6])
            day = int(date[6:8])
            flight_id = text[1]
            sta = text[9]
            dist_m = float(text[4])
            alt = float(text[6])
            x =  float(text[2])  # UTM x-coordinate, meters
            y = float(text[3])  # UTM y-coordinate, meters
            index = stations[stations == sta].index[0]
            elev = float(elevations[index])
            closest_time = float(text[5])
            if (closest_time not in comp_times):
                continue
            index = comp_times.index(closest_time)
            cc,_ = get_speed_of_sound(alt, closest_time, x, y)
            cc_array.append(cc)
            fr_dists.append(abs(np.sqrt(float(text[4])**2 + (float(text[6])-elev)**2)))
            fr_speeds.append(float(text[7]))
            if absolute_time:
                fr_times.append(closest_time + (np.sqrt(float(text[4])**2 + (float(text[6])-elev)**2)/cc) - 1.55e9)
            else:
                tarrive = calc_time(closest_time,dist_m,alt,cc)
                file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/DH8A_data_picks/inversepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight_id) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight_id) + '.csv'

                if Path(file_name).exists():
                    coords = []
                    with open(file_name, 'r') as file:
                        for line in file:
                            pick_data = line.split(',')
                            start_time = float(pick_data[2])
                            break
                fr_times.append(tarrive - start_time)

    scatter1 = axs[ii,0].scatter(inverse_speeds, fr_speeds, c='k', s=15, zorder=2)
    max_min = inverse_speeds + fr_speeds
    axs[ii,0].set_xlim(min(max_min) - 2, max(max_min) + 2)
    axs[ii,0].set_ylim(min(max_min) - 2, max(max_min) + 2)
    if error_bar:
        axs[ii,0].errorbar(inverse_speeds, fr_speeds, xerr=error_vel, fmt='none', c='k', zorder=1)
        max_min_s = inverse_speeds + fr_speeds 
        maxe = max(error_vel)
        axs[ii,0].set_xlim(min(max_min_s) - maxe - 5, max(max_min_s) + maxe + 5)
        axs[ii,0].set_ylim(min(max_min_s) - maxe - 5, max(max_min_s) + maxe + 5)
    axs[0,0].set_title("Velocity (m/s)", fontsize=10)
    axs[ii,0].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[ii,0].set_aspect('equal')
    axs[ii,0].set_xlabel('Inversion Results', fontsize=8)
    axs[ii,0].set_ylabel('Ground Truth', fontsize=8)
    axs[ii,0].tick_params(axis='both', labelsize=8)
    m, b = fit_l1_line(inverse_speeds, fr_speeds)
    x = np.linspace(min(inverse_speeds), max(inverse_speeds), 100)
    axs[ii,0].plot(x, m * x + b, color='k')

    scatter2 = axs[ii,1].scatter(inverse_dists, fr_dists, c='k', s=15, zorder=2)
    max_min = inverse_dists + fr_dists
    axs[ii,1].set_xlim(min(max_min) - 100, max(max_min) + 100)
    axs[ii,1].set_ylim(min(max_min) - 100, max(max_min) + 100)
    if error_bar:
        axs[ii,1].errorbar(inverse_dists, fr_dists, xerr=error_dist, fmt='none', c='k', zorder=1)
        max_min_d = inverse_dists + fr_dists 
        maxe = max(error_dist)
        axs[ii,1].set_xlim(min(max_min_d) - maxe - 100, max(max_min_d) + maxe + 100)
        axs[ii,1].set_ylim(min(max_min_d) - maxe - 100, max(max_min_d) + maxe + 100)
    axs[0,1].set_title("Distance (m)", fontsize=10)
    axs[ii,1].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[ii,1].set_aspect('equal', adjustable='box')
    axs[ii,1].set_xlabel('Inversion Results', fontsize=8)
    axs[ii,1].set_ylabel('Ground Truth', fontsize=8)
    axs[ii,1].tick_params(axis='both', labelsize=8)
    m, b = fit_l1_line(inverse_dists, fr_dists)
    x = np.linspace(min(inverse_dists), max(inverse_dists), 100)
    axs[ii,1].plot(x, m * x + b, color='k')


    scatter3 = axs[ii,2].scatter(c_array, cc_array, c='k', s=15, zorder=2)
    max_min = c_array + cc_array
    axs[ii,2].set_xlim(min(max_min) - 2, max(max_min) + 2)
    axs[ii,2].set_ylim(min(max_min) - 2, max(max_min) + 2)

    if error_bar:
        axs[ii,2].errorbar(c_array, cc_array, xerr=error_time, fmt='none', c='k', zorder=1)
        max_min_c = c_array + cc_array 
        maxe = max(error_c)
        axs[ii,2].set_xlim(min(max_min_c) - maxe - 5, max(max_min_c) + maxe + 5)
        axs[ii,2].set_ylim(min(max_min_c) - maxe - 5, max(max_min_c) + maxe + 5)
    axs[0,2].set_title("Speed of Sound (m/s)", fontsize=10)
    axs[ii,2].set_xlabel('From Inversion', fontsize=8)
    axs[ii,2].set_ylabel('From Modeling', fontsize=8)
    axs[ii,2].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[ii,2].tick_params(axis='both', labelsize=8)
    axs[ii,2].set_aspect('equal', adjustable='box')
    m, b = fit_l1_line(c_array, cc_array)
    x = np.linspace(min(c_array), max(c_array), 100)
    axs[ii,2].plot(x, m * x + b, color='k')

    scatter4 = axs[ii,3].scatter(inverse_times, fr_times, c='k', s=15, zorder=2)
    max_min = inverse_times + fr_times
    if not absolute_time:
        axs[ii,3].set_xlim(min(max_min)-2, max(max_min)+2)
        axs[ii,3].set_ylim(min(max_min)-2, max(max_min)+2)
    if error_bar:
        axs[ii,3].errorbar(inverse_times, fr_times, xerr=error_time, fmt='none', c='k', zorder=1)
        max_min_t = inverse_times + fr_times 
        maxe = max(error_time)
        axs[ii,3].set_xlim(min(max_min_t) - maxe - 20, max(max_min_t) + maxe + 20)
        axs[ii,3].set_ylim(min(max_min_t) - maxe - 20, max(max_min_t) + maxe + 20)

    axs[ii,3].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[0,3].set_title("Time (s)", fontsize=10)
    axs[ii,3].set_xlabel('Inversion Results', fontsize=8)
    axs[ii,3].set_ylabel('Flightradar24 + Modeling', fontsize=8)
    axs[ii,3].tick_params(axis='both', labelsize=8)
    m, b = fit_l1_line(inverse_times, fr_times)
    x = np.linspace(min(inverse_times), max(inverse_times), 100)
    axs[ii,3].plot(x, m * x + b, color='k')
    if absolute_time:
        axs[ii,3].set_xscale('log')
        axs[ii,3].set_yscale('log')
    axs[ii,3].set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.show()
