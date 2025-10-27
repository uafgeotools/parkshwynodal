import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

main_text = False
repo_path = '/home/irseppi/REPOSITORIES/parkshwynodal/'
flightradar_path = '/scratch/irseppi/nodal_data/flightradar24/'
file_in = open(repo_path + 'input/node_crossings_db_UTM.txt','r')

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

col_equip = []
tail = []

tail_num_dict = {}
flight_num_dict = {}
tail_flight_dict = {}

#Get data from crossings database
for text in file_in.readlines():
    lines = text.split(',')
    date = lines[0]
    flight_num = int(lines[1])
    equip = lines[10]
    if equip not in jet + Turboprop + piston + Heli:
        continue
    if equip not in tail_num_dict:
        tail_num_dict[equip] = []
    if equip not in flight_num_dict:
        flight_num_dict[equip] = []
    if flight_num not in flight_num_dict[equip]:
        flight_num_dict[equip].extend([flight_num])

    #Extract tailnumber for flightradar24 data
    flight_data = pd.read_csv(flightradar_path + date + '_flights.csv', sep=",")
    flight = flight_data['flight_id']
    tailnumber = flight_data['aircraft_id']
    for i,fly in enumerate(flight):
        if float(fly) == float(flight_num):
            col_equip.append(equip)
            if int(tailnumber[i]) not in tail_num_dict[equip]:
                tail_num_dict[equip].extend([int(tailnumber[i])])
            if int(tailnumber[i]) not in tail_flight_dict:
                tail_flight_dict[int(tailnumber[i])] = []
            tail_flight_dict[int(tailnumber[i])].extend([int(flight_num)])
            break
file_in.close()

file_jets = open(repo_path + 'output/GT_flight_param_inv_DB.txt', 'r')
file_nonjets = open(repo_path + 'output/NGT_flight_param_inv_DB.txt', 'r')

equip_overtone_dict = {}
equip_count_dict = {}
equip_diff_dict = {}
tail_nums_inverted = {}
line_count_dict = {}

for eq in jet + Turboprop + piston + Heli:
    count1 = 0
    count2 = 0
    equip_overtone_dict[eq] = []
    equip_count_dict[eq] = []
    tail_nums_inverted[eq] = []

    for tt in col_equip:
        if tt == eq:
            count1 += 1
        
    #Get overtone data for jets from the inversion database with the initial model calculated from groundtruth flight parameters
    if eq in jet:
        file_jets = open(repo_path + 'output/GT_flight_param_inv_DB.txt', 'r')
        for line in file_jets:
            lines = line.split(',')
            flight_n = int(lines[1])
            if flight_n not in flight_num_dict[eq] or lines[-2] == "Forward Model":
                continue

            for tail_n, flight_numbers in  tail_flight_dict.items():
                if int(flight_n) in flight_numbers:
                    break
                else:
                    tail_n = None
            if tail_n is None:
                continue

            count2 += 1
            data = np.array(lines[9].strip('[]').split(' '), dtype=float)

            if tail_n not in tail_nums_inverted[eq]:
                tail_nums_inverted[eq].extend([tail_n])
            equip_overtone_dict[eq].extend(data)
        equip_count_dict[eq].extend([count1, count2])
        file_jets.close()

    #Get overtone data for non-jets from the inversion database with the initial model calculated only from the spectrograms
    else:
        file_nonjets = open(repo_path + 'output/NGT_flight_param_inv_DB.txt', 'r')
        for line in file_nonjets:
            lines = line.split(',')
            if lines[0] != eq or lines[-2] == "Forward Model":
                continue
            flight_n = int(lines[2])

            for tail_n, flight_numbers in  tail_flight_dict.items():
                if flight_n in flight_numbers:
                    break
                else:
                    tail_n = None
            if tail_n is None:
                continue

            count2 += 1
            data = np.array(lines[10].strip('[]').split(' '), dtype=float)

            if tail_n not in tail_nums_inverted[eq]:
                tail_nums_inverted[eq].extend([tail_n])
            equip_overtone_dict[eq].extend(data)
        equip_count_dict[eq].extend([count1, count2])
        file_nonjets.close()

    if main_text == True and eq not in jet:
        if eq == 'C185':
            med = 19.5
            line_count = 13
        elif eq == 'C182':
            med = 36
            line_count = 7
        elif eq == 'C206':
            med = 18.7
            line_count = 14
        elif eq == 'DHC2':
            med = 17.5
            line_count = 15
        elif eq == 'GA8':
            med = 20
            line_count = 13
        elif eq == 'PA31':
            med = 18.5
            line_count = 15
        elif eq == 'DH8A':
            med = 60 
            line_count = 16
        elif eq == 'B190':
            med = 24.5
            line_count = 8
        elif eq == 'BE20':
            med = 27.5
            line_count = 8
        elif eq == 'C208':
            med = 85 
            line_count = 9
        elif eq == 'PC12':
            med = 28
            line_count = 8
        elif eq == 'DH3T':
            med = 26
            line_count = 9
        elif eq == 'R44':
            med = 13.35
            line_count = 20
        line_count_dict[eq] = line_count
        equip_diff_dict[eq] = med

    if main_text == False and eq not in jet:
        if eq == 'CH7B':
            med = 19.5
            line_count = 13
        elif eq == 'PA30':
            med = 20.3
            line_count = 13
        elif eq == 'PA32':
            med = 20
            line_count = 14
        elif eq == 'C172':
            med = 19.5
            line_count = 13
        elif eq == 'C180':
            med = 20
            line_count = 12
        elif eq == 'B18T':
            med = 32.3
            line_count = 15
        elif eq == 'C441':
            med = 32
            line_count = 8
        elif eq == 'AT73':
            med = 68.6
            line_count = 15
        elif eq == 'SW4':
            med = 25.7
            line_count = 10
        line_count_dict[eq] = line_count
        equip_diff_dict[eq] = med

title_size = 20
tick_size = 12
text_size = 12

if main_text == True:
    fig, ax = plt.subplots(6, 3, figsize=(20, 24), sharex=True)
    # Track which axes have data
    axes_with_data = set()

    for i, (equip, peaks) in enumerate(equip_overtone_dict.items()):
        equip_count = equip_count_dict[equip]
        if equip not in jet:
            med = equip_diff_dict[equip]
            line_count = line_count_dict[equip]

        if i ==  len(jet) + len(Turboprop):
            label_count = 'crossings: ' + str(equip_count[1]) + '/' + str(equip_count[0])
            label_tail = 'tail numbers: '+ str(len(tail_nums_inverted[equip])) + '/' + str(len(tail_num_dict[equip])) 
        else:
            label_count = str(equip_count[1]) + '/' + str(equip_count[0])
            label_tail = str(len(tail_nums_inverted[equip])) + '/' + str(len(tail_num_dict[equip]))
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
        if equip not in jet:
            med = equip_diff_dict[equip]
            line_count = line_count_dict[equip]
        label_count = str(equip_count[1]) + '/' + str(equip_count[0])
        label_tail = str(len(tail_nums_inverted[equip])) + '/' + str(len(tail_num_dict[equip]))
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

