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
time_in = approach_data.iloc[:, 5]

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
stations = seismo_data['Station']
elevations = seismo_data['Elevation']
file_list = glob.glob('/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results_ngt/*.txt')

fr_dists = []
fr_speeds = []
fr_times = []

inverse_dists = []
inverse_speeds = []
inverse_times = []
for file_name in file_list:
    with open(file_name, 'r') as file:
        equip = file_name.split('/')[-1].split('_')[0]
        if equip in ['B737', 'B738', 'B739', 'B77W', 'B772', 'B788', 'B789', 'B763', 'B744','B733','B732','B77L','B748','CRJ2', 'A332', 'A359', 'E75S']:
            continue
        for line in file.readlines():
            lines = line.split(',')
            date = lines[0]
            month = int(date[4:6])
            day = date[6:8]
            closest_time = float(lines[3])
            comp_time = lines[3]
            flight_num = lines[1]
            sta = lines[2]
            abs_time = float(lines[7])
            if lines[13] == "Forward Model":
                continue
            ins = stations[stations == sta].index[0]
            elev = float(elevations[ins])
            closest_index = None
            file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/inversepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight_num) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight_num) + '.csv'
            if not os.path.exists(file_name):
                continue

            else:
                coords = []
                with open(file_name, 'r') as file:
                    for line in file:
                        pick_data = line.split(',')
                        coords.append((float(pick_data[0]), float(pick_data[1])))
                    if len(pick_data) == 4:
                        start_time = float(pick_data[2])
                    else:
                        file.close() 
                        continue
                file.close()
            if len(coords) == 0:
                continue

            elif equip == 'C185':
                start_time = start_time - 120
            for ii, ss in enumerate(sta_loc):
                if float(ss) == float(sta) and int(flight_id[ii]) == int(flight_num):

                    closest_index = ii
                    fr_dists.append(abs(np.sqrt(float(dist_m[closest_index])**2 + (float(alt[closest_index])-elev)**2)))
                    fr_speeds.append(float(speeds[closest_index]))
                    fr_times.append(closest_time - start_time)

                    inverse_dists.append(abs(float(lines[5])))
                    inverse_speeds.append(abs(float(lines[4])))
                    inverse_times.append(abs_time - start_time)
                    break
            if closest_index is None:
                print(f"Closest time not found for flight {flight_num} at station {sta}")
                continue
inverse_dists = np.array(inverse_dists)
inverse_speeds = np.array(inverse_speeds)
fr_dists = np.array(fr_dists)
fr_speeds = np.array(fr_speeds)

fig, axs = plt.subplots(2, 3, figsize=(10, 10), sharey=False)
axs[0, 0].scatter(inverse_speeds, fr_speeds, c='k', s=15, zorder=2)
axs[0, 0].set_xlim(20,180)
axs[0, 0].set_ylim(20,180)
axs[0, 0].set_aspect('equal')
axs[0, 0].set_title("Speed (m/s)", fontsize=12)
axs[0, 0].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0, 0].set_xlabel('Inversion Results', fontsize=10)
axs[0, 0].set_ylabel('flightradar24', fontsize=10)

#plot text in the top left corner of the first subplot
squared_differences = (np.array(inverse_speeds) - np.array(fr_speeds)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)
rms_speed = rmsd
axs[0, 0].text(0.05, 0.85, 'RMSD = {:.1f}'.format(rmsd), transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')
axs[0, 0].text(0.05, 0.90, '\u03C3 = 30.0', transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')
axs[0, 0].text(0.05, 0.95, 'n = {}'.format(str(np.sum((inverse_speeds > 20) & (inverse_speeds < 180) & (fr_speeds > 20) & (fr_speeds < 180)))), transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')

axs[0, 1].scatter(inverse_dists, fr_dists, c='k', s=15, zorder=2)
axs[0, 1].set_xlim(0, 8500)
axs[0, 1].set_ylim(0, 8500)
axs[0, 1].set_aspect('equal')
axs[0, 1].set_title("Distance (m)", fontsize=12)
axs[0, 1].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0, 1].set_aspect('equal', adjustable='box')
axs[0, 1].set_xlabel('Inversion Results', fontsize=10)
axs[0, 1].set_ylabel('flightradar24', fontsize=10)

squared_differences = (np.array(inverse_dists) - np.array(fr_dists)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)
rms_dist = rmsd
axs[0, 1].text(0.05, 0.85, 'RMSD = {:.1f}'.format(rmsd), transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')
axs[0, 1].text(0.05, 0.90, '\u03C3 = 500.0', transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')
axs[0, 1].text(0.05, 0.95, 'n = {}'.format(str(np.sum((inverse_dists > 0) & (inverse_dists < 8500) & (fr_dists > 0) & (fr_dists < 8500)))), transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')
diff_speed = np.array(inverse_speeds) - np.array(fr_speeds)
diff_dist = np.array(inverse_dists) - np.array(fr_dists)
diff_time = np.array(inverse_times) - np.array(fr_times)

bin = int((np.max(diff_speed) - np.min(diff_speed)) * 4)
mad = np.nanmedian(np.abs(diff_speed - np.median(diff_speed)))
axs[1, 0].hist(diff_speed, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 0].set_xlim(-20,10)
axs[1, 0].set_ylabel(str(np.sum((diff_speed > -20) & (diff_speed < 10))) + '/' + str(len(diff_speed)) + ' samples', fontsize=10)
axs[1, 0].axvline(np.median(diff_speed) - rms_speed, color='red', linestyle='--')
axs[1, 0].axvline(np.median(diff_speed) + rms_speed, color='red', linestyle='--')
axs[1, 0].axvline(np.median(diff_speed), color='red', linestyle='--', linewidth=2)
axs[1, 0].set_title('Median v Difference (m/s): {:.1f} ± {:.1f}'.format(np.median(diff_speed), rms_speed), fontsize=10)
axs[1, 0].set_xlabel('inversion - flightradar24', fontsize=10)

bin = int((np.max(diff_dist) - np.min(diff_dist)) / 27)
mad = np.nanmedian(np.abs(diff_dist - np.median(diff_dist)))
axs[1, 1].set_xlim(-2000,1500)
axs[1, 1].hist(diff_dist, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 1].set_ylabel(str(np.sum((diff_dist > -2000) & (diff_dist < 1500))) + '/' + str(len(diff_dist)) + ' samples', fontsize=10)
axs[1, 1].axvline(np.median(diff_dist) - rms_dist, color='red', linestyle='--')
axs[1, 1].axvline(np.median(diff_dist) + rms_dist, color='red', linestyle='--')
axs[1, 1].axvline(np.median(diff_dist), color='red', linestyle='--', linewidth=2)

axs[1, 1].set_title('Median d\u2080 Difference (m): {:.1f} ± {:.1f}'.format(np.median(diff_dist), rms_dist), fontsize=10)
axs[1, 1].set_xlabel('inversion - flightradar24', fontsize=10)


axs[0, 2].scatter(inverse_times,  fr_times, c='k', s=15, zorder=2)
axs[0, 2].set_xlim(20, 240)
axs[0, 2].set_ylim(20, 240)
axs[0, 2].set_aspect('equal')
axs[0, 2].set_title("Time (s)", fontsize=12)
axs[0, 2].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0, 2].set_aspect('equal', adjustable='box')
axs[0, 2].set_xlabel('Inversion Results', fontsize=10)
axs[0, 2].set_ylabel('flightradar24', fontsize=10)

squared_differences = (np.array(inverse_times) - np.array(fr_times)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)
rms_time = rmsd
axs[0, 2].text(0.05, 0.85, 'RMSD = {:.1f}'.format(rmsd), transform=axs[0, 2].transAxes, fontsize=12, va='top', ha='left')
axs[0, 2].text(0.05, 0.90, '\u03C3 = 30.0', transform=axs[0, 2].transAxes, fontsize=12, va='top', ha='left')
axs[0, 2].text(0.05, 0.95, 'n = {}'.format(str(np.sum((inverse_dists > 0) & (inverse_dists < 8500) & (fr_dists > 0) & (fr_dists < 8500)))), transform=axs[0, 2].transAxes, fontsize=12, va='top', ha='left')

bin = int((np.max(diff_time) - np.min(diff_time)))  # Smaller bin size
mad = np.nanmedian(np.abs(diff_time - np.median(diff_time)))
axs[1, 2].hist(diff_time, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 2].set_xlim(-40,20)
axs[1, 2].axvline(np.median(diff_time) - rms_time, color='red', linestyle='--')
axs[1, 2].axvline(np.median(diff_time) + rms_time, color='red', linestyle='--')
axs[1, 2].axvline(np.median(diff_time), color='red', linestyle='--', linewidth=2)

axs[1, 2].set_title('Median t\u2080 Difference (s): {:.1f} ± {:.1f}'.format(np.median(diff_time), rms_time), fontsize=10)
axs[1, 2].set_xlabel('inversion - flightradar24', fontsize=10)


plt.show()
