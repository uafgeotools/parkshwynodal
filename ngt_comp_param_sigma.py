import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from prelude import *

approach_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', sep=",")
flight_id = approach_data.iloc[:, 1]
x =  approach_data.iloc[:, 2] # UTM x-coordinate, meters
y = approach_data.iloc[:, 3]  # UTM y-coordinate, meters
dist_m = approach_data.iloc[:, 4]
closest_time = approach_data.iloc[:, 5]
alt = approach_data.iloc[:, 6]
speeds = approach_data.iloc[:, 7]
sta_loc = approach_data.iloc[:, 9]

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
stations = seismo_data['Station']
elevations = seismo_data['Elevation']
file_name = 'output/inv_results_no_g_truth/C185_full_inv_results.txt'

error_bar = True
fr_dists = []
fr_speeds = []
cc_array = []

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
color = []
fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

with open(file_name, 'r') as file:
    for line in file.readlines():
        lines = line.split(',')
        comp_time = lines[3]
        flight_num = lines[1]
        sta = lines[2]
        print(lines[13])
        if lines[13] == "Forward Model":
            continue
        ins = stations[stations == sta].index[0]
        elev = float(elevations[ins])
        closest_index = None
        for ii, ss in enumerate(sta_loc):
            if float(ss) == float(sta) and int(flight_id[ii]) == int(flight_num):

                closest_index = ii
                cc,_ = get_speed_of_sound(float(alt[closest_index]), float(closest_time[closest_index]), float(x[closest_index]), float(y[closest_index]))
                cc_array.append(cc)
                fr_dists.append(abs(np.sqrt(float(dist_m[closest_index])**2 + (float(alt[closest_index])-elev)**2)))
                fr_speeds.append(float(speeds[closest_index]))

                inverse_dists.append(abs(float(lines[5])))
                inverse_speeds.append(abs(float(lines[4])))
                c_array.append(float(lines[8]))
                color.append((float(lines[8])-331.3)/0.6)

                if error_bar:
                    error_strs = [e for e in lines[10].strip('[]').split(' ') if e.strip() != '']
                    error = np.array(error_strs) 

                    error_vel.append(float(error[0]) / (len(error) + 4))
                    error_dist.append(float(error[1]) / (len(error) + 4))
                    error_c.append(float(error[3]) / (len(error) + 4))
                break
        if closest_index is None:
            print(f"Closest time not found for flight {flight_num} at station {sta}")
            continue
        
scatter1 = axs[0].scatter(inverse_speeds, fr_speeds, c=color, cmap='coolwarm', s=15, zorder=2)
if error_bar:
    axs[0].errorbar(
        inverse_speeds, fr_speeds, xerr=error_vel,
        fmt='none', ecolor='gray', alpha=0.3, capsize=1, zorder=1, linewidth=0.7
    )
axs[0].set_xlim(40,80)
axs[0].set_ylim(40,80)
axs[0].set_title("Velocity (m/s)", fontsize=10)
axs[0].axline((0, 0), slope=1, color='black', linestyle='--')
axs[0].set_aspect('equal')
axs[0].set_xlabel('Inversion Results', fontsize=8)
axs[0].set_ylabel('flightradar24', fontsize=8)
axs[0].tick_params(axis='both', labelsize=8)
#plot text in the top left corner of the first subplot
squared_differences = (np.array(inverse_speeds) - np.array(fr_speeds)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)
axs[0].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[0].transAxes, fontsize=12, va='top', ha='left')

rms_speed = rmsd
scatter2 = axs[1].scatter(inverse_dists, fr_dists, c=color, cmap='coolwarm', s=15, zorder=2)
if error_bar:
    axs[1].errorbar(inverse_dists, fr_dists, xerr=error_dist, fmt='none', ecolor='gray', alpha=0.3, capsize=1, zorder=1, linewidth=0.7)

axs[1].set_xlim(50, 2500)
axs[1].set_ylim(50, 2500)
axs[1].tick_params(axis='both', labelsize=8)
axs[1].set_title("Distance (m)", fontsize=10)
axs[1].axline((0, 0), slope=1, color='black', linestyle='--')
axs[1].set_aspect('equal', adjustable='box')
axs[1].set_xlabel('Inversion Results', fontsize=8)
axs[1].set_ylabel('flightradar24', fontsize=8)
axs[1].tick_params(axis='both', labelsize=8)
squared_differences = (np.array(inverse_dists) - np.array(fr_dists)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)

rms_dist = rmsd
axs[1].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[1].transAxes, fontsize=12, va='top', ha='left')

scatter3 = axs[2].scatter(c_array, cc_array, c=color, cmap='coolwarm', s=15, zorder=2)
if error_bar:
    axs[2].errorbar(c_array, cc_array, xerr=error_c,  fmt='none', ecolor='gray', alpha=0.3, capsize=1, zorder=1, linewidth=0.7)
#axs[2].set_xlim(290, 335)
#axs[2].set_ylim(290, 335)
axs[2].set_title("Sound Speed(m/s)", fontsize=10)
axs[2].set_xlabel('Inversion Results', fontsize=8)
axs[2].set_ylabel('c(T), T from NCPAG2S', fontsize=8)
axs[2].tick_params(axis='both', labelsize=8)
#axs[2].set_aspect('equal', adjustable='box')
squared_differences = (np.array(c_array) - np.array(cc_array)) ** 2
mean_squared_difference = np.mean(squared_differences)
rmsd = np.sqrt(mean_squared_difference)

rms_c = rmsd

axs[2].axline((0, 0), slope=1, color='black', linestyle='--')
axs[2].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[2].transAxes, fontsize=12, va='top', ha='left')


diff_speed = np.array(inverse_speeds) - np.array(fr_speeds)
diff_dist = np.array(inverse_dists) - np.array(fr_dists)
diff_c = np.array(c_array) - np.array(cc_array)

axs[0].text(0.05, 0.95, '\u03C3 = 30', transform=axs[0].transAxes, fontsize=12, va='top', ha='left')
axs[1].text(0.05, 0.95, '\u03C3 = 500', transform=axs[1].transAxes, fontsize=12, va='top', ha='left')
axs[2].text(0.05, 0.95, '\u03C3 = 100', transform=axs[2].transAxes, fontsize=12, va='top', ha='left')
# Add a single colorbar for the entire figure
cbar = fig.colorbar(scatter3, ax=axs[2], orientation='vertical', pad=0.1)
cbar.set_label('Temperature (°C)')
plt.tight_layout()
plt.show()
plt.close()

fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=False, layout='constrained')
bin = int((np.max(diff_speed) - np.min(diff_speed)) * 3)
axs[0].hist(diff_speed, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[0].set_ylabel(str(len(diff_speed)-1) + '/' + str(len(diff_speed)) + ' samples')
axs[0].axvline(np.mean(diff_speed) - rms_speed, color='red', linestyle='--')
axs[0].axvline(np.mean(diff_speed) + rms_speed, color='red', linestyle='--')
axs[0].axvline(np.mean(diff_speed), color='red', linestyle='--', linewidth=2)
axs[0].set_title('Median Velocity Difference (m/s): {:.2f} ± {:.2f}'.format(np.median(diff_speed), rms_speed), fontsize=14)
axs[0].set_xlabel('inversion - flightradar24')

bin = int((np.max(diff_dist) - np.min(diff_dist)) / 20)
axs[1].hist(diff_dist, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1].set_ylabel(str(len(diff_dist)) + '/' + str(len(diff_dist)) + ' samples')
axs[1].axvline(np.mean(diff_dist) - rms_dist, color='red', linestyle='--')
axs[1].axvline(np.mean(diff_dist) + rms_dist, color='red', linestyle='--')
axs[1].axvline(np.mean(diff_dist), color='red', linestyle='--', linewidth=2)

axs[1].set_title('Median Distance Difference (m): {:.2f} ± {:.2f}'.format(np.median(diff_dist), rms_dist), fontsize=14)
axs[1].set_xlabel('inversion - flightradar24')

bin = int((np.max(diff_c) - np.min(diff_c)) / 3)
axs[2].hist(diff_c, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[2].set_ylabel(str(len(diff_c)-1) + '/' + str(len(diff_c)) + ' samples')
axs[2].axvline(np.mean(diff_c) - rms_c, color='red', linestyle='--')
axs[2].axvline(np.mean(diff_c) + rms_c, color='red', linestyle='--')
axs[2].axvline(np.mean(diff_c), color='red', linestyle='--', linewidth=2)
axs[2].set_title('Median Sound Speed Difference (m/s): {:.2f} ± {:.2f}'.format(np.median(diff_c), rms_c), fontsize=14)
#axs[2].set_yticks(np.arange(1, 9, 1))
axs[2].set_xlabel('inversion - c(T), T from NCPAG2')

plt.show()
