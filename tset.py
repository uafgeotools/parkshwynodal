import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.doppler_funcs import *
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
file_list = glob.glob('/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results_no_g_truth_320/*.txt')


fr_speeds = []
inverse_speeds = []
angles = []

fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharey=False)
for file_name in file_list:
    with open(file_name, 'r') as file:
        equip = file_name.split('/')[-1].split('_')[0]
        if equip in ['B737', 'B738', 'B739', 'B77W', 'B772', 'B788', 'B789', 'B763', 'B744','B733','B732','B77L','B748','CRJ2', 'A332', 'A359', 'E75S']:
            continue
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
                    fr_speeds.append(float(speeds[closest_index]))

                    inverse_speeds.append(abs(float(lines[4])))
                    an = np.arctan(float(dist_m[closest_index]) / (float(alt[closest_index]) - elev))
                    #if -60 > np.rad2deg(an) or np.rad2deg(an) > 60:
                    positive = abs(np.rad2deg(an))
                    print(positive)
                    angles.append(np.cos(np.deg2rad(positive)))
                    #else:
                     #   angles.append(1)

                    break
            if closest_index is None:
                print(f"Closest time not found for flight {flight_num} at station {sta}")
                continue

inverse_speeds = np.array(inverse_speeds)
fr_speeds = np.array(fr_speeds)

axs[0].scatter(inverse_speeds/angles, fr_speeds, c = 'k', s=15, zorder=2)
axs[0].set_title("Velocity (m/s)", fontsize=10)
axs[0].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0].set_aspect('equal')
axs[0].set_xlabel('Inversion Results', fontsize=8)
axs[0].set_ylabel('flightradar24', fontsize=8)
axs[0].tick_params(axis='both', labelsize=8)

squared_differences = (np.array(inverse_speeds) - np.array(fr_speeds)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)

diff_speed = np.array(inverse_speeds) - np.array(fr_speeds)

bin = int((np.max(diff_speed) - np.min(diff_speed)) * 4)
axs[1].hist(diff_speed, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1].set_ylabel(str(np.sum((diff_speed > -20) & (diff_speed < 10))) + '/' + str(len(diff_speed)) + ' samples')
axs[1].axvline(np.mean(diff_speed) - rmsd, color='red', linestyle='--')
axs[1].axvline(np.mean(diff_speed) + rmsd, color='red', linestyle='--')
axs[1].axvline(np.mean(diff_speed), color='red', linestyle='--', linewidth=2)
axs[1].set_title('Median Velocity Difference (m/s): {:.2f} ± {:.2f}'.format(np.median(diff_speed), rmsd), fontsize=10)
axs[1].set_xlabel('inversion - flightradar24')
plt.show()