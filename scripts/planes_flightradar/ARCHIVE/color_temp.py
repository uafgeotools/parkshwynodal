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
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', 'r')
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
file_list = ['C185data_atm_full.txt']

# Define titles and corresponding colors for each word
title = [
    'Full Harmonics/Varying Temp'
]

# Process each file in the list
for gg, file in enumerate(file_list):
    file = open(file, 'r')
    date = []  # List to store indices
    y = 0  # Counter for lines
    time_new = []  # List for new times
    v0_new = []  # List for new velocities
    distance_new = []  # List for new distances
    delf0_new = []  # Placeholder for new frequency differences
    color_temp = []  
    pppp_new = []  # List for new peaks
    med_new = []  # List for median new frequency differences
    mad_new = []  # List for median absolute deviation of new frequency differences
    date_new = []  # List for new peak indices
    date_all = []  # List for all indices
    temp_c = []  # List for temperatures
    error_new = []
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
        time_new.append(float(lines[4]))

        v0_new.append(float(lines[5]))
        flight_num = float(lines[1])
        distance_new.append(float(lines[6]))
        date.append(y)

        peaks_new = np.array(lines[7])
        peaks_new = str(peaks_new)
        peaks_new = np.char.replace(peaks_new, '[', '')
        peaks_new = np.char.replace(peaks_new, ']', '')
        peaks_new = str(peaks_new)
        peaks_new = np.array(peaks_new.split(' '))

        # Calculate median frequency differences for new peaks
        f2 = []
        for i in range(len(peaks_new)):
            peak_new = peaks_new[i]
            pppp_new.append(float(peak_new))
            date_new.append(time)
            color_temp.append(Tc)
            if i == 0:
                continue
            diff = float(peaks_new[i]) - float(peaks_new[i - 1])
            if diff > 22 or diff < 18:
                continue
            f2.append(diff)
        med_new.append(np.nanmedian(f2))
        mad_new.append([np.median(np.absolute(f2 - np.median(f2)))])
        date_all.append(time)

       # Process old and new peaks
        error = np.array(lines[8])
        error = str(error)
        error = np.char.replace(error, '[', '')
        error = np.char.replace(error, ']', '')
        error = str(error)
        error = np.array(error.split(' '))
        for t in range(3,len(error)):
            error_new.append(float(error[t]))
        error_vel = []
        error_dist = []
        error_time = []
        

        error_vel.append(float(error[0]))
        error_dist.append(float(error[1]))
        error_time.append(float(error[2]))

    # Create plots for the current file
    #fig = plt.figure() #figsize=(10, 12))
    #fig.suptitle(title[gg], fontsize=16, y=0.95)
    # Define grid layout for subplots
    #gs = GridSpec(3, 3, figure=fig, width_ratios=[5, 1, 1])

    #ax1 = fig.add_subplot(gs[0, 0])
    #ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    #ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)
    fig,ax1 = plt.subplots(1, 1, sharex=False, figsize=(12,6))     

    ax1.margins(x=0)
    #ax1.grid(axis='both') 
    ax2 = fig.add_axes([0.83, 0.11, 0.07, 0.77], sharey=ax1)
    ax3 = fig.add_axes([0.90, 0.11, 0.07, 0.77], sharey=ax1) 
    # Scatter plot to compare old and new peaks (ax1) and their median frequency differences (ax2)
    ax1.margins(x=0)
    ax1.set_axisbelow(True)
    ax1.scatter(pppp_new, date_new,  c=color_temp,cmap='coolwarm', edgecolors='black', linewidth=0.3, zorder=2)
    ax1.errorbar(pppp_new, date_new, xerr=error_new, fmt='none', c='k', zorder=1)
    ax1.grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    ax1.set_ylabel('Index')
    ax1.set_xlabel('Frequency, Hz ($f_n$)')
    ax1.set_xlim(25, 300)
    ax1.set_xticks(range(0, 275, 1))
    #ax1.text(-0.04, 1.05, 'a)', transform=ax1.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    ax2.set_axisbelow(True)
    ax2.scatter(med_new, date_all, c=temp_c, cmap='coolwarm', edgecolors='black', linewidth=0.3)
    ax2.grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
    ax2.set_xlabel('\u0394' + '$f_{median}$, Hz')
    #ax2.text(-0.1, 1.05, 'b)', transform=ax2.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Add subplot for MAD values
    ax3.set_axisbelow(True)
    ax3.scatter(mad_new, date_all, c=temp_c, cmap='coolwarm', edgecolors='black', linewidth=0.3)
    ax3.grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
    ax3.set_xlabel('\u0394' + '$f_{MAD}$, Hz')
    #ax3.text(-0.1, 1.05, 'c)', transform=ax3.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
    '''
    # Create additional subplots for velocity, distance, and time
    axs = fig.subplots(1, 3, sharey=True) #gridspec_kw={'height_ratios': [1, 1]}, sharey=True)

    # Velocity scatter plots
    axs[0].set_axisbelow(True)
    axs[0].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    scatter = axs[0].scatter((np.array(v0_new) - np.array(v0_old)), date, c=temp_c, s=15, cmap='coolwarm', label='Velocity Residuals')
    axs[0].set_xlabel("$v_0^A - v_0^B$, m/s")
    axs[0].set_xlim(-1.5, 1)
    axs[0].set_ylabel('Index')
    #axs[0].text(-0.1, 1.05, 'g)', transform=axs[1, 0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Distance scatter plots
    axs[1].set_axisbelow(True)
    axs[1].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5)
    scatter = axs[1].scatter((np.array(distance_new) - np.array(distance_old)), date, s=15, c=temp_c, cmap='coolwarm', label='Distance Residuals')
    axs[1].set_xlabel("$l^A - l^B$, m")
    axs[1].set_xlim(-50, 30)
    axs[1].tick_params(left=False, labelleft=False)
    #axs[1].text(-0.1, 1.05, 'h)', transform=axs[1, 1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

    # Time scatter plots
    axs[2].set_axisbelow(True)
    axs[2].grid(which='major', axis='both', color='gray', linestyle='--', linewidth=0.5, zorder=-1)
    scatter = axs[2].scatter((np.array(time_new) - np.array(time_old)), date, s=15, c=temp_c, cmap='coolwarm', label='Time Residuals')
    axs[2].set_xlabel("$t_0^A - t_0^B$, s")
    axs[2].tick_params(left=False, labelleft=False)
    axs[2].set_xlim(-2, 1)
    #axs[2].text(-0.1, 1.05, 'i)', transform=axs[1, 2].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
    '''
    # Add colorbar for temperature
    #fig.colorbar(scatter, ax=axs[1, 2], orientation='vertical', label='Temperature (°C)')

    # Adjust layout and display the plot
    #plt.tight_layout(rect=[0.01, 0.05, 0.95, 0.95], h_pad=0, w_pad=0)
    #plt.savefig(f'/home/irseppi/REPOSITORIES/parkshwynodal/output/{file_list[gg].split(".")[0]}.pdf', dpi=300, bbox_inches='tight')
    plt.show()
