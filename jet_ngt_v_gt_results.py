import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

title_size = 20
tick_size = 12
text_size = 12

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

jet = ['B737', 'B738', 'B739']
fig, ax = plt.subplots(3, 2, figsize=(20, 24), sharex=True)
for i in range(0,2):
    equip_overtone_dict = {}
    equip_count_dict = {}
    equip_diff_dict = {}
    tail_ums_inverted = {}
    line_count_dict = {}

    for eq in jet:
        count1 = 0
        count2 = 0
        if eq not in equip_overtone_dict:
            equip_overtone_dict[eq] = []
            equip_count_dict[eq] = []
            tail_ums_inverted[eq] = []

        for equip in col_equip:
            if equip == eq:
                count1 += 1
        if i == 0:
            # Define the directory where your files are located
            file = '/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results/' + eq + '_full_inv_results.txt' 

            with open(file, 'r') as f:
                # Read the data from the file and append it to the list
                data = []
                fdiff = []
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
                    peak_old = 0
                    for peak in peaks:
                        if peak == '':
                            continue
                        peak = float(peak)
                        if len(peaks) == 0 or peak == peaks[0]:
                            peak_old = float(peak)
                            continue

                        diff = float(peak) - float(peak_old)

                        fdiff.append(diff)
                        peak_old = float(peak)
                        data.append(peak)
                    if tail_n not in tail_ums_inverted[eq]:
                        tail_ums_inverted[eq].append([tail_n])
                equip_overtone_dict[eq].extend(data)
                equip_count_dict[eq].extend([count1, count2])
        else:
            file = '/home/irseppi/REPOSITORIES/parkshwynodal/NGT_flight_param_inv_DB.txt'
            with open(file, 'r') as f:
                # Read the data from the file and append it to the list
                data = []
                fdiff = []
                for line in f.readlines():
                    lines = line.split(',')
                    if lines[0] != eq:
                        continue
                    flight_n = lines[2]
                    for i, flight_num in enumerate(flight_nums):
                        if float(flight_n) == float(flight_num):
                            tail_n = tail[i]
                            break
                    if lines[-2] == "Forward Model" or lines[-5] == "00":
                        continue

                    count2 += 1
                    peaks = np.array(lines[10])
                    peaks = str(peaks) 
                    peaks = peaks.replace('[', '').replace(']', '')
                    peaks = np.array(peaks.split(' '))
                    peak_old = 0
                    for peak in peaks:
                        if peak == '':
                            continue
                        peak = float(peak)
                        if len(peaks) == 0 or peak == peaks[0]:
                            peak_old = float(peak)
                            continue

                        diff = float(peak) - float(peak_old)

                        fdiff.append(diff)
                        peak_old = float(peak)
                        data.append(peak)
                    if tail_n not in tail_ums_inverted[eq]:
                        tail_ums_inverted[eq].append([tail_n])
                equip_overtone_dict[eq].extend(data)
                equip_count_dict[eq].extend([count1, count2])
    
        med = np.median(fdiff)
        if eq == 'B739':
            med = 35.5
            line_count = 8

        elif eq == 'B738':
            med = 73
            line_count = 3

        elif eq == 'B737':
            med = 68
            line_count = 3

        line_count_dict[eq] = line_count

        equip_diff_dict[eq] = med





    for j, (equip, peaks) in enumerate(equip_overtone_dict.items()):
        print(j, equip, len(peaks))
        equip_count = equip_count_dict[equip]
        med = equip_diff_dict[equip]
        line_count = line_count_dict[equip]


        label_count = str(equip_count[1]) + '/' + str(equip_count[0])
        label_tail = str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip]))

        if j == 0 and i == 0:
            ax[j, i].set_title('Inversion results with ground truth', fontsize=title_size, fontweight='bold')
        elif j == 0 and i == 1:
            ax[j, i].set_title('Inversion results without ground truth', fontsize=title_size, fontweight='bold')
        bins = np.arange(min(peaks), max(peaks) + 3, 3)
        ax[j, i].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
        ax[j, i].text(0.99, 0.95, equip , transform=ax[j, i].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        ax[j, i].text(0.99, 0.85, label_count, transform=ax[j, i].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        ax[j, i].text(0.99, 0.75, label_tail, transform=ax[j, i].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        ax[j, i].text(0.99, 0.65, len(peaks), transform=ax[j, i].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        counts, _ = np.histogram(peaks, bins=bins)
        ax[j, i].set_yticks([0,counts.max()])
        ax[j, i].tick_params(axis='y', labelsize=tick_size)

# Set x-ticks for specific axes and ensure they are visible
for (row, col) in [(2,1), (2,0)]:
    ax[row, col].set_xticks(np.arange(0, 300, 25))
    ax[row, col].tick_params(axis='x', which='both', length=3, labelbottom=True, labelsize=tick_size)
    ax[row, col].set_xlabel('Frequency (Hz)', fontsize=text_size)

plt.xlim(5, 300)
plt.tight_layout(pad=2.7, w_pad=0.5, h_pad=0)


plt.show()

