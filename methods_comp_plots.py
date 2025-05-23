import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj
from matplotlib.gridspec import GridSpec
from matplotlib import colors as mcolors


# Initialize UTM projection for coordinate conversion
utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

# Read input file containing station crossing data
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt', 'r')
latc = []  # List to store latitudes
lonc = []  # List to store longitudes
timec = []  # List to store times
altc = []  # List to store altitudes

# Parse the input file line by line
for line in file_in.readlines():
    text = line.split(',')

    # Extract time and UTM coordinates
    timec.append(float(text[5]))
    x = float(text[2])
    y = float(text[3])

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)
    latc.append(lat)
    lonc.append(lon)
    altc.append(float(text[4]) * 0.0003048)  # Convert altitude from feet to kilometers
file_in.close()

# List of files to process and their corresponding titles
file_list = ['1o_atmc_v_2c.txt', '2c_1o_v_full.txt', 'full_atmc_v_2c.txt', 'atmc_1o_v_full.txt']

# Define titles and corresponding colors for each word
title = [
    'One Harmonic/Fixed Temp vs. One Harmonic/Varying Temp',
    'One Harmonic/Fixed Temp vs. Full Harmonics/Fixed Temp',
    'Full Harmonics/Fixed Temp vs. Full Harmonics/Varying Temp',
    'One Harmonic/Varying Temp vs. Full Harmonics/Varying Temp'
]

# Process each file in the list
for gg, file in enumerate(file_list):
    file = open(file, 'r')
    date = []  # List to store indices
    y = 0  # Counter for lines
    time_new = []  # List for new times
    time_old = []  # List for old times

    v0_old = []  # List for old velocities
    v0_new = []  # List for new velocities

    distance_old = []  # List for old distances
    distance_new = []  # List for new distances

    delf0_new = []  # Placeholder for new frequency differences
    delf0_old = []  # Placeholder for old frequency differences

    pppp_old = []  # List for old peaks
    pppp_new = []  # List for new peaks
    med_old = []  # List for median old frequency differences
    med_new = []  # List for median new frequency differences
    mad_old = []  # List for median absolute deviation of old frequency differences
    mad_new = []  # List for median absolute deviation of new frequency differences
    date_old = []  # List for old peak indices
    date_new = []  # List for new peak indices
    date_all = []  # List for all indices
    temp_c = []  # List for temperatures

    # Parse the file line by line
    for line in file.readlines():
        y += 1
        counts = []

        lines = line.split(',')
        time = float(lines[3])

        # Match the time with the corresponding latitude, longitude, and altitude
        for i in range(len(latc)):
            if timec[i] == time:
                lat = latc[i]
                lon = lonc[i]
                alt = altc[i]
                break

        # Construct the input file path for atmospheric data
        input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'
        try:
            file = open(input_files, 'r')
        except:
            continue
        data = json.load(file)

        # Extract metadata and data from the atmosphere JSON file
        metadata = data['metadata']
        sourcefile = metadata['sourcefile']
        datetim = metadata['time']['datetime']
        latitude = metadata['location']['latitude']
        longitude = metadata['location']['longitude']
        parameters = metadata['parameters']

        data_list = data['data']

        # Convert data to a pandas DataFrame
        data_frame = pd.DataFrame(data_list)

        # Find the closest altitude index
        z_index = None
        hold = np.inf
        for item in data_list:
            if item['parameter'] == 'Z':
                for i in range(len(item['values'])):
                    if abs(float(item['values'][i]) - float(alt)) < hold:
                        hold = abs(float(item['values'][i]) - float(alt))
                        z_index = i

        # Extract temperature at the closest altitude
        for item in data_list:
            if item['parameter'] == 'T':
                Tc = -273.15 + float(item['values'][z_index])  # Convert from Kelvin to Celsius
        temp_c.append(Tc)

        # Extract other parameters from the line
        time_old.append(float(lines[4]))
        time_new.append(float(lines[9]))
        v0_old.append(float(lines[5]))
        v0_new.append(float(lines[10]))
        flight_num = float(lines[1])
        distance_old.append(float(lines[6]))
        distance_new.append(float(lines[11]))
        date.append(y)

        # Process old and new peaks
        peaks_old = np.array(lines[7])
        peaks_old = str(peaks_old)
        peaks_old = np.char.replace(peaks_old, '[', '')
        peaks_old = np.char.replace(peaks_old, ']', '')
        peaks_old = str(peaks_old)
        peaks_old = np.array(peaks_old.split(' '))

        peaks_new = np.array(lines[12])
        peaks_new = str(peaks_new)
        peaks_new = np.char.replace(peaks_new, '[', '')
        peaks_new = np.char.replace(peaks_new, ']', '')
        peaks_new = str(peaks_new)
        peaks_new = np.array(peaks_new.split(' '))

        # Calculate median frequency differences for old peaks
        f1 = []
        for i in range(len(peaks_old)):
            peak_old = peaks_old[i]
            pppp_old.append(float(peak_old))
            date_old.append(y)
            if i == 0:
                continue
            diff = float(peaks_old[i]) - float(peaks_old[i - 1])
            if diff > 22 or diff < 18:
                continue
            f1.append(diff)
        med_old.append(np.nanmedian(f1))
        mad_old.append([np.median(np.absolute(f1 - np.median(f1)))])
        # Calculate median frequency differences for new peaks
        f2 = []
        for i in range(len(peaks_new)):
            peak_new = peaks_new[i]
            pppp_new.append(float(peak_new))
            date_new.append(y)
            if i == 0:
                continue
            diff = float(peaks_new[i]) - float(peaks_new[i - 1])
            if diff > 22 or diff < 18:
                continue
            f2.append(diff)
        med_new.append(np.nanmedian(f2))
        mad_new.append([np.median(np.absolute(f2 - np.median(f2)))])
        date_all.append(y)

    # Create plots for the current file
    fig = plt.figure(figsize=(10, 12))
    fig.suptitle(title[gg], fontsize=16, y=0.95)
    plt.figtext(0.40, 0.91, "(A: Square)", fontsize=16, color='cyan', weight="bold", ha='right')
    plt.figtext(0.60, 0.91, "(B: Triangle)", fontsize=16, color='orange', weight="bold", ha='left')
    # Define grid layout for subplots
    gs = GridSpec(3, 3, figure=fig, width_ratios=[5, 1, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)
    # Scatter plot to compare old and new peaks (ax1) and their median frequency differences (ax2)
    ax1.margins(x=0)
    ax1.set_axisbelow(True)
    ax1.scatter(pppp_old, date_old, c='cyan', s=15, marker='s', edgecolors='black', linewidth=0.3)
    ax1.scatter(pppp_new, date_new, c='orange', s=15, marker='^', edgecolors='black', linewidth=0.3)
    ax1.grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    ax1.set_ylabel('Index')
    ax1.set_xlabel('Frequency, Hz ($f_n$)')
    ax1.set_xlim(25, 275)
    ax1.set_xticks(range(0, 251, 25))
    ax1.text(-0.04, 1.05, 'a)', transform=ax1.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    ax2.set_axisbelow(True)
    ax2.scatter(med_old, date_all, c='cyan', s=15, marker='s', edgecolors='black', linewidth=0.3)
    ax2.scatter(med_new, date_all, c='orange', s=15, marker='^', edgecolors='black', linewidth=0.3)
    ax2.grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
    ax2.set_xlabel('\u0394' + '$f_{median}$, Hz')
    ax2.text(-0.1, 1.05, 'b)', transform=ax2.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Add subplot for MAD values
    ax3.set_axisbelow(True)
    ax3.scatter(mad_old, date_all, c='cyan', s=15, marker='s', edgecolors='black', linewidth=0.3)
    ax3.scatter(mad_new, date_all, c='orange', s=15, marker='^', edgecolors='black', linewidth=0.3)
    ax3.grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
    ax3.set_xlabel('\u0394' + '$f_{MAD}$, Hz')
    ax3.text(-0.1, 1.05, 'c)', transform=ax3.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Create additional subplots for velocity, distance, and time
    axs = fig.subplots(2, 3, gridspec_kw={'height_ratios': [1, 1], 'top': 0.59, 'hspace': 0.3, 'wspace': 0.15})

    # Velocity scatter plots
    axs[0, 0].set_axisbelow(True)
    axs[0, 0].scatter(v0_old, date, c='cyan', s=15, marker='s', edgecolors='black', linewidth=0.3)
    axs[0, 0].scatter(v0_new, date, c='orange', s=15, marker='^', edgecolors='black', linewidth=0.3)
    axs[0, 0].set_xlabel('Velocity, m/s ($v_0$)')
    axs[0, 0].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    axs[0, 0].set_xlim(45, 85)
    axs[0, 0].set_ylabel('Index')
    axs[0, 0].text(-0.1, 1.05, 'd)', transform=axs[0, 0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    axs[1, 0].set_axisbelow(True)
    axs[1, 0].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    scatter = axs[1, 0].scatter((np.array(v0_new) - np.array(v0_old)), date, c=temp_c, s=15, cmap='coolwarm', label='Velocity Residuals')
    axs[1, 0].set_xlabel("$v_0^A - v_0^B$, m/s")
    axs[1, 0].set_xlim(-1.5, 1)
    axs[1, 0].set_ylabel('Index')
    axs[1, 0].text(-0.1, 1.05, 'g)', transform=axs[1, 0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Distance scatter plots
    axs[0, 1].set_axisbelow(True)
    axs[0, 1].scatter(distance_old, date, c='cyan', s=15, marker='s', edgecolors='black', linewidth=0.3)
    axs[0, 1].scatter(distance_new, date, c='orange', s=15, marker='^', edgecolors='black', linewidth=0.3)
    axs[0, 1].set_xlabel('Distance, m ($l$)')
    axs[0, 1].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    axs[0, 1].set_xlim(0, 3500)
    axs[0, 1].tick_params(left=False, labelleft=False)
    axs[0, 1].text(-0.1, 1.05, 'e)', transform=axs[0, 1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    axs[1, 1].set_axisbelow(True)
    axs[1, 1].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    scatter = axs[1, 1].scatter((np.array(distance_new) - np.array(distance_old)), date, s=15, c=temp_c, cmap='coolwarm', label='Distance Residuals')
    axs[1, 1].set_xlabel("$l^A - l^B$, m")
    axs[1, 1].set_xlim(-50, 30)
    axs[1, 1].tick_params(left=False, labelleft=False)
    axs[1, 1].text(-0.1, 1.05, 'h)', transform=axs[1, 1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Time scatter plots
    axs[0, 2].set_axisbelow(True)
    axs[0, 2].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    axs[0, 2].set_xlim(110, 120)
    axs[0, 2].scatter(time_old, date, c='cyan', s=15, marker='s', edgecolors='black', linewidth=0.3)
    axs[0, 2].scatter(time_new, date, c='orange', s=15, marker='^', edgecolors='black', linewidth=0.3)
    axs[0, 2].set_xlabel('Time, s ($t_0$)')
    axs[0, 2].tick_params(left=False, labelleft=False)
    axs[0, 2].text(-0.1, 1.05, 'f)', transform=axs[0, 2].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    axs[1, 2].set_axisbelow(True)
    axs[1, 2].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5, zorder=-1)
    scatter = axs[1, 2].scatter((np.array(time_new) - np.array(time_old)), date, s=15, c=temp_c, cmap='coolwarm', label='Time Residuals')
    axs[1, 2].set_xlabel("$t_0^A - t_0^B$, s")
    axs[1, 2].tick_params(left=False, labelleft=False)
    axs[1, 2].set_xlim(-2, 1)
    axs[1, 2].text(-0.1, 1.05, 'i)', transform=axs[1, 2].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Add colorbar for temperature
    fig.colorbar(scatter, ax=axs[1, 2], orientation='vertical', label='Temperature (°C)')

    # Adjust layout and display the plot
    plt.tight_layout(rect=[0.01, 0.05, 0.95, 0.95], h_pad=0, w_pad=0)
    plt.savefig(f'/home/irseppi/REPOSITORIES/parkshwynodal/output/{file_list[gg].split(".")[0]}.pdf', dpi=300, bbox_inches='tight')
    #plt.show()
