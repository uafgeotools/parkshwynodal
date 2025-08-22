import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from prelude import *
import glob

approach_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', sep=",")
flight_id = approach_data.iloc[:, 1]
dist_m = approach_data.iloc[:, 4]
alt = approach_data.iloc[:, 6]
speeds = approach_data.iloc[:, 7]
sta_loc = approach_data.iloc[:, 9]

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
stations = seismo_data['Station']
elevations = seismo_data['Elevation']
file_list = glob.glob('output/inv_results_no_g_truth/*.txt')
print(file_list)
error_bar = False
fr_dists = []
fr_speeds = []

inverse_dists = []
inverse_speeds = []

if error_bar:
    error_vel = []
    error_dist = []

fig, axs = plt.subplots(2, 2, figsize=(10, 10), sharey=False)
for file_name in file_list:
    with open(file_name, 'r') as file:
        for line in file.readlines():
            lines = line.split(',')
            comp_time = lines[3]
            flight_num = lines[1]
            sta = lines[2]
            if lines[13] == "Forward Model":
                continue
            ins = stations[stations == sta].index[0]
            elev = float(elevations[ins])
            closest_index = None
            for ii, ss in enumerate(sta_loc):
                if float(ss) == float(sta) and int(flight_id[ii]) == int(flight_num):

                    closest_index = ii
                    fr_dists.append(abs(np.sqrt(float(dist_m[closest_index])**2 + (float(alt[closest_index])-elev)**2)))
                    fr_speeds.append(float(speeds[closest_index]))

                    inverse_dists.append(abs(float(lines[5])))
                    inverse_speeds.append(abs(float(lines[4])))

                    if error_bar:
                        error_strs = [e for e in lines[10].strip('[]').split(' ') if e.strip() != '']
                        error = np.array(error_strs) 

                        error_vel.append(float(error[0]) / (len(error) + 4))
                        error_dist.append(float(error[1]) / (len(error) + 4))
                    break
            if closest_index is None:
                print(f"Closest time not found for flight {flight_num} at station {sta}")
                continue
inverse_dists = np.array(inverse_dists)
inverse_speeds = np.array(inverse_speeds)
fr_dists = np.array(fr_dists)
fr_speeds = np.array(fr_speeds)
axs[0, 0].scatter(inverse_speeds, fr_speeds, c='k', s=15, zorder=2)
if error_bar:
    axs[0].errorbar(
        inverse_speeds, fr_speeds, xerr=error_vel,
        fmt='none', ecolor='gray', alpha=0.3, capsize=1, zorder=1, linewidth=0.7
    )
axs[0, 0].set_xlim(20,180)
axs[0, 0].set_ylim(20,180)
axs[0, 0].set_title("Velocity (m/s)", fontsize=10)
axs[0, 0].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0, 0].set_aspect('equal')
axs[0, 0].set_xlabel('Inversion Results', fontsize=8)
axs[0, 0].set_ylabel('flightradar24', fontsize=8)
axs[0, 0].tick_params(axis='both', labelsize=8)
#plot text in the top left corner of the first subplot
squared_differences = (np.array(inverse_speeds) - np.array(fr_speeds)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)
axs[0, 0].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')

rms_speed = rmsd
axs[0, 1].scatter(inverse_dists, fr_dists, c='k', s=15, zorder=2)
if error_bar:
    axs[0, 1].errorbar(inverse_dists, fr_dists, xerr=error_dist, fmt='none', ecolor='gray', alpha=0.3, capsize=1, zorder=1, linewidth=0.7)

axs[0, 1].set_xlim(0, 8500)
axs[0, 1].set_ylim(0, 8500)
axs[0, 1].tick_params(axis='both', labelsize=8)
axs[0, 1].set_title("Distance (m)", fontsize=10)
axs[0, 1].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0, 1].set_aspect('equal', adjustable='box')
axs[0, 1].set_xlabel('Inversion Results', fontsize=8)
axs[0, 1].set_ylabel('flightradar24', fontsize=8)
axs[0, 1].tick_params(axis='both', labelsize=8)
squared_differences = (np.array(inverse_dists) - np.array(fr_dists)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)

rms_dist = rmsd
axs[0, 1].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')

diff_speed = np.array(inverse_speeds) - np.array(fr_speeds)
diff_dist = np.array(inverse_dists) - np.array(fr_dists)

axs[0, 0].text(0.05, 0.95, '\u03C3 = 30', transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')
axs[0, 1].text(0.05, 0.95, '\u03C3 = 500', transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')

bin = int((np.max(diff_speed) - np.min(diff_speed)) * 4)
axs[1, 0].hist(diff_speed, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 0].set_xlim(-20,10)
axs[1, 0].set_ylabel(str(np.sum((diff_speed > -20) & (diff_speed < 10))) + '/' + str(len(diff_speed)) + ' samples')
axs[1, 0].axvline(np.mean(diff_speed) - rms_speed, color='red', linestyle='--')
axs[1, 0].axvline(np.mean(diff_speed) + rms_speed, color='red', linestyle='--')
axs[1, 0].axvline(np.mean(diff_speed), color='red', linestyle='--', linewidth=2)
axs[1, 0].set_title('Median Velocity Difference (m/s): {:.2f} ± {:.2f}'.format(np.median(diff_speed), rms_speed), fontsize=10)
axs[1, 0].set_xlabel('inversion - flightradar24')

bin = int((np.max(diff_dist) - np.min(diff_dist)) / 27)
axs[1, 1].set_xlim(-2000,1500)
axs[1, 1].hist(diff_dist, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 1].set_ylabel(str(np.sum((diff_dist > -2000) & (diff_dist < 1500))) + '/' + str(len(diff_dist)) + ' samples')
axs[1, 1].axvline(np.mean(diff_dist) - rms_dist, color='red', linestyle='--')
axs[1, 1].axvline(np.mean(diff_dist) + rms_dist, color='red', linestyle='--')
axs[1, 1].axvline(np.mean(diff_dist), color='red', linestyle='--', linewidth=2)

axs[1, 1].set_title('Median Distance Difference (m): {:.2f} ± {:.2f}'.format(np.median(diff_dist), rms_dist), fontsize=10)
axs[1, 1].set_xlabel('inversion - flightradar24')

plt.show()
