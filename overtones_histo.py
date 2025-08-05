import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
col_equip = []
flight_nums = []
tail = []
tail_num_dict = {}
for text in file_in.readlines():
    lines = text.split(',')
    equip = lines[10]
    date = lines[0]
    flight_num = lines[1]
    col_equip.append(equip)
    flight_nums.append(flight_num)
    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + date + '_flights.csv', sep=",")
    flight = flight_data['flight_id']
    tailnumber = flight_data['aircraft_id']
    if equip not in tail_num_dict:
        tail_num_dict[equip] = []
    for i,fly in enumerate(flight):
        if float(fly) == float(flight_num):
            tail.append(tailnumber[i])
            if tailnumber[i] not in tail_num_dict[equip]:
                tail_num_dict[equip].append(tailnumber[i])
            break
file_in.close()
jet = ['B737', 'B738', 'B739', 'B77W']
Turboprop = ['DH8A','B190','BE20','C208','PC12','DH3T']
piston = ['C185','C182','C206','DHC2','GA8','PA31']
Heli = ['R44']
# Move to suplumentary material
# piston = ['CH7B', 'PA30', 'PA32', 'C172','C180']
# Turboprop = ['B18T','C441','AT73','SW4']
# jet = ['B733', 'B763', 'B772', 'B788', 'B789']
equip_overtone_dict = {}
equip_count_dict = {}
tail_ums_inverted = {}
for eq in jet + Turboprop + piston + Heli:
    count1 = 0
    count2 = 0
    if eq not in equip_overtone_dict:
        equip_overtone_dict[eq] = []
        equip_count_dict[eq] = []
        tail_ums_inverted[eq] = []

    for equip in col_equip:
        if equip == eq:
            count1 += 1
    # Define the directory where your files are located
    file = 'output/inv_results/' + eq + '_full_inv_results.csv' 

    with open(file, 'r') as f:
        # Read the data from the file and append it to the list
        data = []
        for line in f.readlines():
            lines = line.split(',')
            flight_n = lines[1]
            for i, flight_num in enumerate(flight_nums):
                if float(flight_n) == float(flight_num):
                    tail_n = tail[i]
                    break
            if lines[-2] == "Forward Model" or lines[-5] == "00":
                continue

            count2 += 1
            peaks = np.array(lines[9])
            peaks = str(peaks) 
            peaks = peaks.replace('[', '').replace(']', '')
            peaks = np.array(peaks.split(' '))
            for peak in peaks:
                if peak == '':
                    continue
                peak = float(peak)
            
                data.append(peak)
            if tail_n not in tail_ums_inverted[eq]:
                tail_ums_inverted[eq].append([tail_n])
        equip_overtone_dict[eq].extend(data)
        equip_count_dict[eq].extend([count1, count2])
fig, ax = plt.subplots(6, 3, figsize=(20, 25), sharex=True)

# Track which axes have data
axes_with_data = set()

for i, (equip, peaks) in enumerate(equip_overtone_dict.items()):
    equip_count = equip_count_dict[equip]
    label_count = 'Crossings: ' + str(equip_count[1]) + '/' + str(equip_count[0])
    label_tail = 'Tailnumbers: '+ str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip])) 
    if equip in jet:
        if i == 0:
            ax[i, 2].set_title('Jet Aircrafts', fontsize=14, fontweight='bold')
        bins = np.arange(min(peaks), max(peaks) + 3, 3)
        ax[i, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
        ax[i, 2].text(0.99, 0.95, equip, transform=ax[i, 2].transAxes, fontsize=10, va='top', ha='right')
        ax[i, 2].text(0.99, 0.85, label_count, transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right')
        ax[i, 2].text(0.99, 0.75, label_tail, transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right')
        ax[i, 2].text(0.99, 0.65, len(peaks), transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right')
        counts, _ = np.histogram(peaks, bins=bins)
        ax[i, 2].set_yticks([0,counts.max()])
        axes_with_data.add((i, 2))
    elif equip in Turboprop:
        idx = i - len(jet)
        if i == 9:
            ax[0, 1].set_title('Turboprop Aircrafts', fontsize=14, fontweight='bold')
        bins = np.arange(min(peaks), max(peaks) + 3, 3)
        ax[idx, 1].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
        ax[idx, 1].text(0.99, 0.95, equip, transform=ax[idx, 1].transAxes, fontsize=10, va='top', ha='right')
        ax[idx, 1].text(0.99, 0.85, label_count, transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right')
        ax[idx, 1].text(0.99, 0.75, label_tail, transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right')
        ax[idx, 1].text(0.99, 0.65, len(peaks), transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right')
        counts, _ = np.histogram(peaks, bins=bins)
        ax[idx, 1].set_yticks([0,counts.max()])
        axes_with_data.add((idx, 1))
    elif equip in piston:
        idx = i - len(jet) - len(Turboprop)
        if i == len(jet) + len(Turboprop):
            ax[idx, 0].set_title('Piston Aircrafts', fontsize=14, fontweight='bold')
        bins = np.arange(min(peaks), max(peaks) + 3, 3)
        ax[idx, 0].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
        ax[idx, 0].text(0.99, 0.95, equip, transform=ax[idx, 0].transAxes, fontsize=10, va='top', ha='right')
        ax[idx, 0].text(0.99, 0.85, label_count, transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right')
        ax[idx, 0].text(0.99, 0.75, label_tail, transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right')
        ax[idx, 0].text(0.99, 0.65, len(peaks), transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right')
        counts, _ = np.histogram(peaks, bins=bins)
        ax[idx, 0].set_yticks([0,counts.max()])
        axes_with_data.add((idx, 0))
    elif equip in Heli:
        ax[-1, 2].set_title('Helicopter (Piston)', fontsize=14, fontweight='bold')
        bins = np.arange(min(peaks), max(peaks) + 3, 3)
        ax[-1, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
        ax[-1, 2].text(0.99, 0.95, equip, transform=ax[-1, 2].transAxes, fontsize=10, va='top', ha='right')
        ax[-1, 2].text(0.99, 0.85, label_count, transform=ax[-1, 2].transAxes, fontsize=9, va='top', ha='right')
        ax[-1, 2].text(0.99, 0.75, label_tail, transform=ax[-1, 2].transAxes, fontsize=9, va='top', ha='right')
        ax[-1, 2].text(0.99, 0.65, len(peaks), transform=ax[-1, 2].transAxes, fontsize=9, va='top', ha='right')
        counts, _ = np.histogram(peaks, bins=bins)
        ax[-1, 2].set_yticks([0,counts.max()])
        axes_with_data.add((5, 2))  


# Remove outline for axes with no data
for row in range(ax.shape[0]):
    for col in range(ax.shape[1]):
        if (row, col) not in axes_with_data:
            for name, spine in ax[row, col].spines.items():
                #if name != 'top':
                spine.set_visible(False)
            ax[row, col].set_xticks([])
            ax[row, col].set_yticks([])
            ax[row, col].set_facecolor('none')  # Remove plot background
            plt.setp(ax[row, col].get_xticklabels(), visible=False)
            ax[row, col].tick_params(axis='x', which='both', length=0, labelbottom=False)
        else:
            plt.setp(ax[row, col].get_yticklabels(), visible=False)
            ax[row, col].tick_params(axis='y', which='both', length=3, labelleft=False)

# Set x-ticks for specific axes and ensure they are visible
for (row, col) in [(5,1), (5,0), (3,2), (5,2)]:
    ax[row, col].set_xticks(np.arange(0, 300, 25))
    ax[row, col].tick_params(axis='x', which='both', length = 3, labelbottom=True)

#plt.subplots_adjust(hspace=0, wspace=0.1)  # Small vertical and horizontal space between subplots
plt.xlim(5, 300)
plt.tight_layout(pad=2, w_pad=0.5, h_pad=0)

fig.savefig('histogram.png', dpi=300, bbox_inches='tight')
plt.show()

