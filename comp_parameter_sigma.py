import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.doppler_funcs import *

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
stations = seismo_data['Station']
elevations = seismo_data['Elevation']
file_names = ['FIXED_C_50_100_1000_30_1_DH8A.txt','50_100_1000_30_1_DH8A.txt','50_10_200_30_80_DH8A.txt','50_1_1_30_160_DH8A.txt']

fr_dists = []
fr_speeds = []
cc_array = []
error_bar = False

fig, axs = plt.subplots(4, 3, figsize=(10, 15), sharey=False, layout='constrained')
for ii, ff in enumerate(file_names):
    file = open(ff,'r')
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
            cc,_ = get_speed_of_sound(alt, closest_time, x, y)
            cc_array.append(cc)
            fr_dists.append(abs(np.sqrt(float(text[4])**2 + (float(text[6])-elev)**2)))
            fr_speeds.append(float(text[7]))

    scatter1 = axs[ii,0].scatter(inverse_speeds, fr_speeds, c='k', s=15, zorder=2)
    axs[ii,0].set_xlim(90,160)
    axs[ii,0].set_ylim(90,160)
    axs[ii,0].set_xticks(np.arange(100, 160, 10))
    axs[ii,0].set_yticks(np.arange(100, 160, 10))
    axs[0,0].set_title("Velocity (m/s)", fontsize=10)
    axs[ii,0].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[ii,0].set_aspect('equal')
    axs[ii,0].set_xlabel('Inversion Results', fontsize=8)
    axs[ii,0].set_ylabel('flightradar24', fontsize=8)
    axs[ii,0].tick_params(axis='both', labelsize=8)
    #plot text in the top left corner of the first subplot
    squared_differences = (np.array(inverse_speeds) - np.array(fr_speeds)) ** 2
    mean_squared_difference = np.mean(squared_differences)
    rmsd = np.sqrt(mean_squared_difference)
    axs[ii,0].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[ii,0].transAxes, fontsize=12, va='top', ha='left')
    if ii == 2:
        rms_speed = rmsd
    scatter2 = axs[ii,1].scatter(inverse_dists, fr_dists, c='k', s=15, zorder=2)
    axs[ii,1].set_xlim(4500, 7500)
    axs[ii,1].set_ylim(4500, 7500)
    axs[ii,1].set_xticks(np.arange(5000, 7500, 500))
    axs[ii,1].set_yticks(np.arange(5000, 7500, 500))
    axs[ii,1].tick_params(axis='both', labelsize=8)
    axs[0,1].set_title("Distance (m)", fontsize=10)
    axs[ii,1].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[ii,1].set_aspect('equal', adjustable='box')
    axs[ii,1].set_xlabel('Inversion Results', fontsize=8)
    axs[ii,1].set_ylabel('flightradar24', fontsize=8)
    axs[ii,1].tick_params(axis='both', labelsize=8)
    squared_differences = (np.array(inverse_dists) - np.array(fr_dists)) ** 2
    mean_squared_difference = np.mean(squared_differences)
    rmsd = np.sqrt(mean_squared_difference)
    if ii == 2:
        rms_dist = rmsd
    axs[ii,1].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[ii,1].transAxes, fontsize=12, va='top', ha='left')

    scatter3 = axs[ii,2].scatter(c_array, cc_array, c='k', s=15, zorder=2)
    axs[ii,2].set_xlim(245, 340)
    axs[ii,2].set_ylim(245, 340)
    axs[ii,2].set_xticks(np.arange(250, 340, 10))
    axs[ii,2].set_yticks(np.arange(250, 340, 10))
    axs[0,2].set_title("Sound Speed(m/s)", fontsize=10)
    axs[ii,2].set_xlabel('Inversion Results', fontsize=8)
    axs[ii,2].set_ylabel('c(T), T from NCPAG2S', fontsize=8)
    axs[ii,2].tick_params(axis='both', labelsize=8)
    axs[ii,2].set_aspect('equal', adjustable='box')
    squared_differences = (np.array(c_array) - np.array(cc_array)) ** 2
    mean_squared_difference = np.mean(squared_differences)
    rmsd = np.sqrt(mean_squared_difference)
    if ii == 2:
        rms_c = rmsd
    if ii != 0:
        axs[ii,2].axline((0, 0), slope=1, color='black', linestyle='--')
        axs[ii,2].text(0.05, 0.85, 'RMSD = {:.2f}'.format(rmsd), transform=axs[ii,2].transAxes, fontsize=12, va='top', ha='left')
    else:
        axs[ii,2].axvline(310.72, color='black', linestyle='--')
        axs[ii,2].text(0.05, 0.85, 'RMSD = 0.00', transform=axs[ii,2].transAxes, fontsize=12, va='top', ha='left')
    if ii == 2:
        diff_speed = np.array(inverse_speeds) - np.array(fr_speeds)
        diff_dist = np.array(inverse_dists) - np.array(fr_dists)
        diff_c = np.array(c_array) - np.array(cc_array)
    if ii == 0 or ii == 1:
        axs[ii,0].text(0.05, 0.95, '\u03C3 = 100', transform=axs[ii,0].transAxes, fontsize=12, va='top', ha='left')
        axs[ii,1].text(0.05, 0.95, '\u03C3 = 1000', transform=axs[ii,1].transAxes, fontsize=12, va='top', ha='left')
        axs[ii,2].text(0.05, 0.95, '\u03C3 = 1', transform=axs[ii,2].transAxes, fontsize=12, va='top', ha='left')
    elif ii == 2:
        axs[ii,0].text(0.05, 0.95, '\u03C3 = 10', transform=axs[ii,0].transAxes, fontsize=12, va='top', ha='left')
        axs[ii,1].text(0.05, 0.95, '\u03C3 = 200', transform=axs[ii,1].transAxes, fontsize=12, va='top', ha='left')
        axs[ii,2].text(0.05, 0.95, '\u03C3 = 80', transform=axs[ii,2].transAxes, fontsize=12, va='top', ha='left')
    elif ii == 3:
        axs[ii,0].text(0.05, 0.95, '\u03C3 = 1', transform=axs[ii,0].transAxes, fontsize=12, va='top', ha='left')
        axs[ii,1].text(0.05, 0.95, '\u03C3 = 1', transform=axs[ii,1].transAxes, fontsize=12, va='top', ha='left')
        axs[ii,2].text(0.05, 0.95, '\u03C3 = 160', transform=axs[ii,2].transAxes, fontsize=12, va='top', ha='left')

plt.tight_layout()
fig.savefig("compare_sigma.pdf", dpi=500)
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
axs[0].set_xlim(-1.5, 2.5)
axs[0].set_ylim(0,8.5)
axs[0].set_yticks(np.arange(1, 9, 1))
axs[0].set_xticks(np.arange(-1, 3, 1))
axs[0].set_xlabel('inversion - flightradar24')
bin = int((np.max(diff_dist) - np.min(diff_dist)) / 3)
axs[1].hist(diff_dist, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[1].set_ylabel(str(len(diff_dist)) + '/' + str(len(diff_dist)) + ' samples')
axs[1].axvline(np.mean(diff_dist) - rms_dist, color='red', linestyle='--')
axs[1].axvline(np.mean(diff_dist) + rms_dist, color='red', linestyle='--')
axs[1].axvline(np.mean(diff_dist), color='red', linestyle='--', linewidth=2)
axs[1].set_ylim(0,8.5)
axs[1].set_yticks(np.arange(1, 9, 1))
axs[1].set_title('Median Distance Difference (m): {:.2f} ± {:.2f}'.format(np.median(diff_dist), rms_dist), fontsize=14)
axs[1].set_xlabel('inversion - flightradar24')

bin = int((np.max(diff_c) - np.min(diff_c)) / 3)
axs[2].hist(diff_c, bins=bin, color='k', edgecolor='black', alpha=0.5)
axs[2].set_ylabel(str(len(diff_c)-1) + '/' + str(len(diff_c)) + ' samples')
axs[2].axvline(np.mean(diff_c) - rms_c, color='red', linestyle='--')
axs[2].axvline(np.mean(diff_c) + rms_c, color='red', linestyle='--')
axs[2].axvline(np.mean(diff_c), color='red', linestyle='--', linewidth=2)
axs[2].set_xlim(-12, 22)
axs[2].set_ylim(0,8.5)
axs[2].set_title('Median Sound Speed Difference (m/s): {:.2f} ± {:.2f}'.format(np.median(diff_c), rms_c), fontsize=14)
axs[2].set_yticks(np.arange(1, 9, 1))
axs[2].set_xlabel('inversion - c(T), T from NCPAG2')
fig.savefig("sigma_comp_histo.pdf", dpi=500)
plt.show()