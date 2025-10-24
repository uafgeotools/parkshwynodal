import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

repo_path = '/home/irseppi/REPOSITORIES/parkshwynodal/'
approach_data = pd.read_csv(repo_path + 'input/node_crossings_db_UTM.txt', sep=",")
flight_id = approach_data.iloc[:, 1]
dist_m = approach_data.iloc[:, 4]
alt = approach_data.iloc[:, 6]
speeds = approach_data.iloc[:, 7]
sta_loc = approach_data.iloc[:, 9]
time_in = approach_data.iloc[:, 5]

seismo_data = pd.read_csv(repo_path + 'input/parkshwy_nodes.txt', sep="|")
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

inversion_results = pd.read_csv(repo_path + 'output/NGT_flight_param_inv_DB.txt', sep=",")
equipment = inversion_results['Equipment']
closest_times = inversion_results['Closest_Approach_Timestamp']
flight_numbers = inversion_results['Flight_Number']
sensors = inversion_results['Station']
meas_absolute_times = inversion_results['Meas_Closest_Approach_Timestamp']
data_misfit = inversion_results['Data_Misfit_Value']
meas_speeds = inversion_results['Meas_Speed']
meas_dists = inversion_results['Meas_Distance']

fr_dists = []
fr_speeds = []
fr_times = []

inverse_dists = []
inverse_speeds = []
inverse_times = []

for idx, equip in enumerate(equipment):
    if equip in ['B737', 'B738', 'B739', 'B77W', 'B772', 'B788', 'B789', 'B763', 'B744','B733','B732','B77L','B748','CRJ2', 'A332', 'A359', 'E75S']:
        continue
    closest_time = float(closest_times[idx])
    flight_num = flight_numbers[idx]
    sta = sensors[idx]
    abs_time = float(meas_absolute_times[idx])
    if data_misfit[idx] == "Forward Model":
        continue
    ins = stations[stations == sta].index[0]
    elev = float(elevations[ins])
    closest_index = None
    for ii, ss in enumerate(sta_loc):
        if float(ss) == float(sta) and int(flight_id[ii]) == int(flight_num):
            closest_index = ii
            fr_dists.append(abs(np.sqrt(float(dist_m[closest_index])**2 + (float(alt[closest_index])-elev)**2)))
            fr_speeds.append(float(speeds[closest_index]))
            fr_times.append(closest_time) 

            inverse_dists.append(abs(float(meas_dists[idx])))
            inverse_speeds.append(abs(float(meas_speeds[idx])))
            inverse_times.append(abs_time)
            break
    if closest_index is None:
        print(f"Closest time not found for flight {flight_num} at station {sta}")
        continue
inverse_dists = np.array(inverse_dists)
inverse_speeds = np.array(inverse_speeds)
inverse_times = np.array(inverse_times)
fr_dists = np.array(fr_dists)
fr_speeds = np.array(fr_speeds)
fr_times = np.array(fr_times)

diff_speed = inverse_speeds - fr_speeds
diff_dist = inverse_dists - fr_dists
diff_time = inverse_times - fr_times

mad_speed = np.nanmedian(np.abs(diff_speed - np.median(diff_speed)))
mad_dist = np.nanmedian(np.abs(diff_dist - np.median(diff_dist)))
mad_time = np.nanmedian(np.abs(diff_time - np.median(diff_time)))

fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharey=False)

axs[0, 0].scatter(inverse_speeds, fr_speeds, c='k', s=15, zorder=2)
axs[0, 0].set_xlim(20,180)
axs[0, 0].set_ylim(20,180)
axs[0, 0].set_aspect('equal')
axs[0, 0].set_title("Aircraft Speed (m/s)", fontsize=12)
axs[0, 0].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0, 0].set_xlabel('inversion results', fontsize=10)
axs[0, 0].set_ylabel('flightradar24', fontsize=10)
axs[0, 0].text(0.05, 0.75, 'MAD = {:.1f}'.format(mad_speed), transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')
axs[0, 0].text(0.05, 0.85, '\u03C3 = 30.0', transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')
axs[0, 0].text(0.05, 0.95, 'n = {}'.format(str(np.sum((inverse_speeds > 20) & (inverse_speeds < 240) & (fr_speeds > 20) & (fr_speeds < 240)))), transform=axs[0, 0].transAxes, fontsize=12, va='top', ha='left')

axs[0, 1].scatter(inverse_dists, fr_dists, c='k', s=15, zorder=2)
axs[0, 1].set_xlim(0, 8500)
axs[0, 1].set_ylim(0, 8500)
axs[0, 1].set_aspect('equal')
axs[0, 1].set_title("Distance (m)", fontsize=12)
axs[0, 1].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0, 1].set_aspect('equal', adjustable='box')
axs[0, 1].set_xlabel('inversion results', fontsize=10)
axs[0, 1].set_ylabel('flightradar24', fontsize=10)
axs[0, 1].text(0.05, 0.75, 'MAD = {:.1f}'.format(mad_dist), transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')
axs[0, 1].text(0.05, 0.85, '\u03C3 = 500.0', transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')
axs[0, 1].text(0.05, 0.95, 'n = {}'.format(str(np.sum((inverse_dists > 0) & (inverse_dists < 8500) & (fr_dists > 0) & (fr_dists < 8500)))), transform=axs[0, 1].transAxes, fontsize=12, va='top', ha='left')

bin = int((np.max(diff_speed) - np.min(diff_speed)) * 4)
axs[1, 0].hist(diff_speed, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 0].set_xlim(-20,10)
axs[1, 0].set_ylabel(str(np.sum((diff_speed > -20) & (diff_speed < 10))) + '/' + str(len(diff_speed)) + ' samples', fontsize=10)
axs[1, 0].axvline(np.median(diff_speed) - mad_speed, color='red', linestyle='--')
axs[1, 0].axvline(np.median(diff_speed) + mad_speed, color='red', linestyle='--')
axs[1, 0].axvline(np.median(diff_speed), color='red', linestyle='--', linewidth=2)
axs[1, 0].set_title('Median v: {:.1f} ± {:.1f} m/s'.format(np.median(diff_speed), mad_speed), fontsize=10)
axs[1, 0].set_xlabel('inversion - flightradar24', fontsize=10)

bin = int((np.max(diff_dist) - np.min(diff_dist)) / 27)
axs[1, 1].set_xlim(-2000,1500)
axs[1, 1].hist(diff_dist, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 1].set_ylabel(str(np.sum((diff_dist > -2000) & (diff_dist < 1500))) + '/' + str(len(diff_dist)) + ' samples', fontsize=10)
axs[1, 1].axvline(np.median(diff_dist) - mad_dist, color='red', linestyle='--')
axs[1, 1].axvline(np.median(diff_dist) + mad_dist, color='red', linestyle='--')
axs[1, 1].axvline(np.median(diff_dist), color='red', linestyle='--', linewidth=2)
axs[1, 1].set_title('Median d\u2080: {:.1f} ± {:.1f} m'.format(np.median(diff_dist), mad_dist), fontsize=10)
axs[1, 1].set_xlabel('inversion - flightradar24', fontsize=10)

bin = int((np.max(diff_time) - np.min(diff_time)))  # Smaller bin size
axs[0, 2].hist(diff_time, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[0, 2].set_xlim(-40,20)
axs[0, 2].axvline(np.median(diff_time) - mad_time, color='red', linestyle='--')
axs[0, 2].axvline(np.median(diff_time) + mad_time, color='red', linestyle='--')
axs[0, 2].axvline(np.median(diff_time), color='red', linestyle='--', linewidth=2)
axs[0, 2].text(0.05, 0.95, '\u03C3 = 30.0', transform=axs[0, 2].transAxes, fontsize=12, va='top', ha='left')
axs[0, 2].set_ylabel(str(np.sum((diff_time > -40) & (diff_time < 20))) + '/' + str(len(diff_time)) + ' samples', fontsize=10)
axs[0, 2].set_title('Median t\u2080: {:.1f} ± {:.1f} s'.format(np.median(diff_time), mad_time), fontsize=12)
axs[0, 2].set_xlabel('inversion - flightradar24', fontsize=10)

bin = int((np.max(diff_time) - np.min(diff_time)))  # Smaller bin size
axs[1, 2].hist(diff_time, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1, 2].set_xlim(-40,20)
axs[1, 2].axvline(np.median(diff_time) - mad_time, color='red', linestyle='--')
axs[1, 2].axvline(np.median(diff_time) + mad_time, color='red', linestyle='--')
axs[1, 2].axvline(np.median(diff_time), color='red', linestyle='--', linewidth=2)
axs[1, 2].set_ylabel(str(np.sum((diff_time > -40) & (diff_time < 20))) + '/' + str(len(diff_time)) + ' samples', fontsize=10)
axs[1, 2].set_title('Median t\u2080: {:.1f} ± {:.1f} s'.format(np.median(diff_time), mad_time), fontsize=10)
axs[1, 2].set_xlabel('inversion - flightradar24', fontsize=10)

plt.show()
