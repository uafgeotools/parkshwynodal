import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

main_text = True
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
if main_text == True:
    jet = ['B737', 'B738', 'B739']
    Turboprop = ['DH8A','B190','BE20','C208','PC12','DH3T']
    piston = ['C185','C182','C206','DHC2','GA8','PA31']
    Heli = ['R44']
else:
    # Move to suplumentary material
    piston = ['CH7B', 'PA30', 'PA32', 'C172','C180']
    Turboprop = ['B18T','C441','AT73','SW4']
    jet = ['B733', 'B763', 'B772', 'B77W', 'B788', 'B789']
    Heli = []
equip_overtone_dict = {}
equip_count_dict = {}
equip_diff_dict = {}
tail_ums_inverted = {}
line_count_dict = {}
bc_dict = {}
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
    if eq in jet:
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
    if main_text == True:
        if eq == 'C185':
            med = 19.5
            line_count = 13
            blade_count = '2/3'
        elif eq == 'C182':
            med = 36
            line_count = 7
            blade_count = '2/3'
        elif eq == 'C206':
            med = 18.7
            line_count = 14
            blade_count = '2/3/5'
        elif eq == 'DHC2':
            med = 17.5
            line_count = 15
            blade_count = '2/3'
        elif eq == 'GA8':
            med = 20
            line_count = 13
            blade_count = '2/3'
        elif eq == 'PA31':
            med = 18.5
            line_count = 15
            blade_count = '2/3/4'
        elif eq == 'DH8A':
            med = 60 
            line_count = 16
            blade_count = '4'
        elif eq == 'B190':
            med = 24.5
            line_count = 8
            blade_count = '4/5'
        elif eq == 'BE20':
            med = 27.5
            line_count = 8
            blade_count = '3/4/5'
        elif eq == 'C208':
            med = 85 
            line_count = 9
            blade_count = '3/4/5'
        elif eq == 'PC12':
            med = 28
            line_count = 8
            blade_count = '4/5/7'
        elif eq == 'DH3T':
            med = 26
            line_count = 9
            blade_count = '4'
        elif eq == 'R44':
            med = 13.35
            line_count = 20
            blade_count = '2'
        elif eq == 'B739':
            med = 35.5
            line_count = 8
            blade_count = '24'
        elif eq == 'B738':
            med = 73
            line_count = 3
            blade_count = '24'
        elif eq == 'B737':
            med = 68
            line_count = 3
            blade_count = '18/24'
        line_count_dict[eq] = line_count
        bc_dict[eq] = blade_count
        equip_diff_dict[eq] = med
    if main_text == False:
        if eq == 'CH7B':
            med = 19.5
            line_count = 13
            blade_count = '2'
        elif eq == 'PA30':
            med = 20.3
            line_count = 13
            blade_count = '2/3'
        elif eq == 'PA32':
            med = 20
            line_count = 14
            blade_count = '2/3'
        elif eq == 'C172':
            med = 19.5
            line_count = 13
            blade_count = '2/3'
        elif eq == 'C180':
            med = 20
            line_count = 12
            blade_count = '2/3'
        elif eq == 'B18T':
            med = 32.3
            line_count = 15
            blade_count = '3'
        elif eq == 'C441':
            med = 32
            line_count = 8
            blade_count = '3/4/5'
        elif eq == 'AT73':
            med = 68.6
            line_count = 15
            blade_count = '6'
        elif eq == 'SW4':
            med = 25.7
            line_count = 10
            blade_count = '3/4/5'
        elif eq == 'B733':
            med = 76
            line_count = 3
            blade_count = '18/24/36'
        elif eq == 'B763':
            med = 35
            line_count = 3
            blade_count = '33/38'
        elif eq == 'B772':
            med = 30
            line_count = 9
            blade_count = '22/26'
        elif eq == 'B77W':
            med = 13.35
            line_count = 20
            blade_count = '22'
        elif eq == 'B788':
            med = 35.5
            line_count = 8
            blade_count = '18'
        elif eq == 'B789':
            med = 35
            line_count = 3
            blade_count = '18/22'

        line_count_dict[eq] = line_count
        bc_dict[eq] = blade_count
        equip_diff_dict[eq] = med
if main_text == True:
    fig, ax = plt.subplots(6, 3, figsize=(20, 24), sharex=True)

    # Track which axes have data
    axes_with_data = set()

    for i, (equip, peaks) in enumerate(equip_overtone_dict.items()):
        equip_count = equip_count_dict[equip]
        med = equip_diff_dict[equip]
        line_count = line_count_dict[equip]
        blade_count = bc_dict[equip]
        if i ==  len(jet) + len(Turboprop):
            label_count = 'crossings: ' + str(equip_count[1]) + '/' + str(equip_count[0])
            label_tail = 'tail numbers: '+ str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip])) 
        else:
            label_count = str(equip_count[1]) + '/' + str(equip_count[0])
            label_tail = str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip]))
        if equip in jet:
            if i == 0:
                ax[i, 2].set_title('Jet Aircraft', fontsize=title_size, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[i, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[i, 2].text(0.99, 0.95, equip , transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[i, 2].text(0.99, 0.85, label_count, transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[i, 2].text(0.99, 0.75, label_tail, transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[i, 2].text(0.99, 0.65, len(peaks), transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            counts, _ = np.histogram(peaks, bins=bins)
            ax[i, 2].set_yticks([0,counts.max()])
            ax[i, 2].tick_params(axis='y', labelsize=tick_size)
            axes_with_data.add((i, 2))
        elif equip in Turboprop:
            idx = i - len(jet)
            if i == 8:
                ax[0, 1].set_title('Turboprop Aircraft', fontsize=title_size, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 1].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 1].text(0.99, 0.95, equip , transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.85, label_count, transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.75, label_tail, transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.65, len(peaks), transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            # Show one decimal if med is not an integer, else show ".0"
            med_display = f"{med:.1f}" if float(med).is_integer() else str(round(med, 1))
            ax[idx, 1].text(
                0.99, 0.55, med_display,
                transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
            )
            for g in range(0,line_count):
                ax[idx, 1].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
                if equip == 'DH8A':
                    ax[idx, 1].axvline(x= (1 + g) * 15, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1, alpha=0.3)
                elif equip == 'C208':
                    ax[idx, 1].axvline(x= (1 + g) * 28.333, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1, alpha=0.3)    
            counts, _ = np.histogram(peaks, bins=bins)
            ax[idx, 1].set_yticks([0,counts.max()])
            ax[idx, 1].tick_params(axis='y', labelsize=tick_size)
            axes_with_data.add((idx, 1))
        elif equip in piston:
            idx = i - len(jet) - len(Turboprop)
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 0].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 0].text(0.99, 0.95, equip , transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[idx, 0].text(0.99, 0.85, label_count, transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[idx, 0].text(0.99, 0.75, label_tail, transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            if i == len(jet) + len(Turboprop):
                ax[idx, 0].set_title('Piston Aircraft', fontsize=title_size, fontweight='bold')
                ax[idx, 0].text(0.99, 0.65, 'f\u209B count: ' + str(len(peaks)), transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
                ax[idx, 0].text(0.99, 0.55, '\u0394f\u209B: ' + str(round(med,1)) + ' Hz', transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            else:
                ax[idx, 0].text(0.99, 0.65, len(peaks), transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
                med_display = f"{med:.1f}" if float(med).is_integer() else str(round(med, 1))
                ax[idx, 0].text(0.99, 0.55, med_display, transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            for g in range(0, line_count):
                ax[idx, 0].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
            counts, _ = np.histogram(peaks, bins=bins)
            ax[idx, 0].set_yticks([0,counts.max()])
            ax[idx, 0].tick_params(axis='y', labelsize=tick_size)
            axes_with_data.add((idx, 0))
        elif equip in Heli:
            ax[-2, 2].set_title('Helicopter (Piston)', fontsize=title_size, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[-2, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[-2, 2].text(0.99, 0.95, equip , transform=ax[-2, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[-2, 2].text(0.99, 0.85, label_count, transform=ax[-2, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[-2, 2].text(0.99, 0.75, label_tail, transform=ax[-2, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[-2, 2].text(0.99, 0.65, len(peaks), transform=ax[-2, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax[-2, 2].text(0.99, 0.55, round(med,1), transform=ax[-2, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

            for g in range(0,20):
                ax[-2, 2].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
            counts, _ = np.histogram(peaks, bins=bins)
            ax[-2, 2].set_yticks([0,counts.max()])
            ax[-2, 2].tick_params(axis='y', labelsize=tick_size)
            axes_with_data.add((4, 2))


    # Remove outline for axes with no data
    for row in range(ax.shape[0]):
        for col in range(ax.shape[1]):
            if (row, col) not in axes_with_data:
                for name, spine in ax[row, col].spines.items():
                    spine.set_visible(False)
                ax[row, col].set_xticks([])
                ax[row, col].set_yticks([])
                ax[row, col].set_facecolor('none')  # Remove plot background
                plt.setp(ax[row, col].get_xticklabels(), visible=False)
                ax[row, col].tick_params(axis='x', which='both', length=0, labelbottom=False)

    # Set x-ticks for specific axes and ensure they are visible
    for (row, col) in [(5,1), (5,0), (2,2), (4,2)]:
        ax[row, col].set_xticks(np.arange(0, 300, 25))
        ax[row, col].tick_params(axis='x', which='both', length=3, labelbottom=True, labelsize=tick_size)
        ax[row, col].set_xlabel('Frequency (Hz)', fontsize=text_size)

    plt.xlim(5, 300)
    plt.tight_layout(pad=2.7, w_pad=0.5, h_pad=0)


if main_text == False:

    fig, ax = plt.subplots(6, 3, figsize=(20, 24), sharex=True)

    # Track which axes have data
    axes_with_data = set()

    for i, (equip, peaks) in enumerate(equip_overtone_dict.items()):
        equip_count = equip_count_dict[equip]
        med = equip_diff_dict[equip]
        line_count = line_count_dict[equip]
        blade_count = bc_dict[equip]
        label_count = str(equip_count[1]) + '/' + str(equip_count[0])
        label_tail = str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip]))
        if equip in jet:
            if i == 0:
                ax[i, 2].set_title('Jet Aircraft', fontsize=title_size, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[i, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[i, 2].text(0.99, 0.95, equip , transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[i, 2].text(0.99, 0.85, label_count, transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[i, 2].text(0.99, 0.75, label_tail, transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[i, 2].text(0.99, 0.65, str(len(peaks)), transform=ax[i, 2].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            counts, _ = np.histogram(peaks, bins=bins)
            ax[i, 2].set_yticks([0,counts.max()])
            axes_with_data.add((i, 2))
        elif equip in Turboprop:
            idx = i - len(jet)
            if i == len(jet):
                ax[0, 1].set_title('Turboprop Aircraft', fontsize=title_size, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 1].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 1].text(0.99, 0.95, equip , transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.85, label_count, transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.75, label_tail, transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.65, str(len(peaks)), transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            for g in range(0,line_count):
                ax[idx, 1].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
                if equip == 'AT73':
                    ax[idx, 1].axvline(x= (1 + g) * 17, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1, alpha=0.3)
            med_display = f"{med:.1f}" if float(med).is_integer() else str(round(med, 1))
            ax[idx, 1].text(
                0.99, 0.55, med_display,
                transform=ax[idx, 1].transAxes, fontsize=text_size, va='top', ha='right',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
            )
            counts, _ = np.histogram(peaks, bins=bins)
            ax[idx, 1].set_yticks([0,counts.max()])
            axes_with_data.add((idx, 1))
        elif equip in piston:
            idx = i - len(jet) - len(Turboprop)
            if i == len(jet) + len(Turboprop):
                ax[idx, 0].set_title('Piston Aircraft', fontsize=title_size, fontweight='bold')


            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 0].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 0].text(0.99, 0.95, equip , transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 0].text(0.99, 0.85, label_count, transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 0].text(0.99, 0.75, label_tail, transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 0].text(
                0.99, 0.65, str(len(peaks)),
                transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right',
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
            )
            for g in range(0,line_count):
                ax[idx, 0].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
            med_display = f"{med:.1f}" if float(med).is_integer() else str(round(med, 1))
            ax[idx, 0].text(
                0.99, 0.55, med_display,
                transform=ax[idx, 0].transAxes, fontsize=text_size, va='top', ha='right',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
            )
            counts, _ = np.histogram(peaks, bins=bins)
            ax[idx, 0].set_yticks([0,counts.max()])
            axes_with_data.add((idx, 0))

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

    # Set x-ticks for specific axes and ensure they are visible
    for (row, col) in [(3,1), (4,0), (5,2)]:
        ax[row, col].set_xticks(np.arange(0, 300, 25))
        ax[row, col].tick_params(axis='x', which='both', length = 3, labelbottom=True)
        ax[row, col].set_xlabel('Frequency (Hz)', fontsize=text_size)

    plt.xlim(5, 300)
    plt.tight_layout(pad=2.7, w_pad=0.5, h_pad=0)

plt.show()

