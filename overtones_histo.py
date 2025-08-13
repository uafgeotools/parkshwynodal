import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
main_text = True
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
    med = np.median(fdiff)
    equip_diff_dict[eq] = med
if main_text == True:
    fig, ax = plt.subplots(6, 3, figsize=(20, 25), sharex=True)

    # Track which axes have data
    axes_with_data = set()

    for i, (equip, peaks) in enumerate(equip_overtone_dict.items()):
        equip_count = equip_count_dict[equip]
        med = equip_diff_dict[equip]
        if i ==  len(jet) + len(Turboprop):
            label_count = 'crossings: ' + str(equip_count[1]) + '/' + str(equip_count[0])
            label_tail = 'tail numbers: '+ str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip])) 
        else:
            label_count = str(equip_count[1]) + '/' + str(equip_count[0])
            label_tail = str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip]))
        if equip in jet:
            if i == 0:
                ax[i, 2].set_title('Jet Aircrafts', fontsize=14, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[i, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[i, 2].text(0.99, 0.95, equip, transform=ax[i, 2].transAxes, fontsize=10, va='top', ha='right')
            ax[i, 2].text(0.99, 0.85, label_count, transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right')
            ax[i, 2].text(0.99, 0.75, label_tail, transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right')
            ax[i, 2].text(0.99, 0.65, len(peaks), transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right')
            for g in range(0,20):
                ax[i, 2].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
            counts, _ = np.histogram(peaks, bins=bins)
            ax[i, 2].set_yticks([0,counts.max()])
            axes_with_data.add((i, 2))
        elif equip in Turboprop:
            idx = i - len(jet)
            if i == 8:
                ax[0, 1].set_title('Turboprop Aircrafts', fontsize=14, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 1].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 1].text(0.99, 0.95, equip, transform=ax[idx, 1].transAxes, fontsize=10, va='top', ha='right')
            ax[idx, 1].text(0.99, 0.85, label_count, transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right')
            ax[idx, 1].text(0.99, 0.75, label_tail, transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right')
            ax[idx, 1].text(0.99, 0.65, len(peaks), transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right')
            for g in range(0,20):
                ax[idx, 1].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
            counts, _ = np.histogram(peaks, bins=bins)
            ax[idx, 1].set_yticks([0,counts.max()])
            axes_with_data.add((idx, 1))
        elif equip in piston:
            idx = i - len(jet) - len(Turboprop)
            if i == len(jet) + len(Turboprop):
                ax[idx, 0].set_title('Piston Aircrafts', fontsize=14, fontweight='bold')
                ax[idx, 0].text(0.99, 0.65, 'f\u2080 count: ' + str(len(peaks)), transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right')
            else:
                ax[idx, 0].text(0.99, 0.65, len(peaks), transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 0].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 0].text(0.99, 0.95, equip, transform=ax[idx, 0].transAxes, fontsize=10, va='top', ha='right')
            ax[idx, 0].text(0.99, 0.85, label_count, transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right')
            ax[idx, 0].text(0.99, 0.75, label_tail, transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right')
            for g in range(0,20):
                ax[idx, 0].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
            counts, _ = np.histogram(peaks, bins=bins)
            ax[idx, 0].set_yticks([0,counts.max()])
            axes_with_data.add((idx, 0))
        elif equip in Heli:
            ax[-2, 2].set_title('Helicopter (Piston)', fontsize=14, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[-2, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[-2, 2].text(0.99, 0.95, equip, transform=ax[-2, 2].transAxes, fontsize=10, va='top', ha='right')
            ax[-2, 2].text(0.99, 0.85, label_count, transform=ax[-2, 2].transAxes, fontsize=9, va='top', ha='right')
            ax[-2, 2].text(0.99, 0.75, label_tail, transform=ax[-2, 2].transAxes, fontsize=9, va='top', ha='right')
            ax[-2, 2].text(0.99, 0.65, len(peaks), transform=ax[-2, 2].transAxes, fontsize=9, va='top', ha='right')
            for g in range(0,20):
                ax[-2, 2].axvline(x= (1 + g) * med, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
            counts, _ = np.histogram(peaks, bins=bins)
            ax[-2, 2].set_yticks([0,counts.max()])
            axes_with_data.add((4, 2))  


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
    for (row, col) in [(5,1), (5,0), (2,2), (4,2)]:
        ax[row, col].set_xticks(np.arange(0, 300, 25))
        ax[row, col].tick_params(axis='x', which='both', length = 3, labelbottom=True)

    #plt.subplots_adjust(hspace=0, wspace=0.1)  # Small vertical and horizontal space between subplots
    plt.xlim(5, 300)
    plt.tight_layout(pad=2, w_pad=0.5, h_pad=0)

    fig.savefig('histogram.png', dpi=300, bbox_inches='tight')
    plt.show()
if main_text == False:
    fig, ax = plt.subplots(6, 3, figsize=(17, 25), sharex=True)

    # Track which axes have data
    axes_with_data = set()

    for i, (equip, peaks) in enumerate(equip_overtone_dict.items()):
        equip_count = equip_count_dict[equip]
        label_count = str(equip_count[1]) + '/' + str(equip_count[0])
        label_tail = str(len(tail_ums_inverted[equip])) + '/' + str(len(tail_num_dict[equip]))
        if equip in jet:
            if i == 0:
                ax[i, 2].set_title('Jet Aircrafts', fontsize=14, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[i, 2].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[i, 2].text(0.99, 0.95, equip, transform=ax[i, 2].transAxes, fontsize=10, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[i, 2].text(0.99, 0.85, label_count, transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[i, 2].text(0.99, 0.75, label_tail, transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

            ax[i, 2].text(0.99, 0.65, str(len(peaks)), transform=ax[i, 2].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            counts, _ = np.histogram(peaks, bins=bins)
            ax[i, 2].set_yticks([0,counts.max()])
            axes_with_data.add((i, 2))
        elif equip in Turboprop:
            idx = i - len(jet)
            if i == len(jet):
                ax[0, 1].set_title('Turboprop Aircrafts', fontsize=14, fontweight='bold')
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 1].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 1].text(0.99, 0.95, equip, transform=ax[idx, 1].transAxes, fontsize=10, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.85, label_count, transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.75, label_tail, transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 1].text(0.99, 0.65, str(len(peaks)), transform=ax[idx, 1].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            counts, _ = np.histogram(peaks, bins=bins)
            ax[idx, 1].set_yticks([0,counts.max()])
            axes_with_data.add((idx, 1))
        elif equip in piston:
            idx = i - len(jet) - len(Turboprop)
            if i == len(jet) + len(Turboprop):
                ax[idx, 0].set_title('Piston Aircrafts', fontsize=14, fontweight='bold')

            ax[idx, 0].text(
                0.99, 0.65, str(len(peaks)),
                transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right',
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
            )
            bins = np.arange(min(peaks), max(peaks) + 3, 3)
            ax[idx, 0].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
            ax[idx, 0].text(0.99, 0.95, equip, transform=ax[idx, 0].transAxes, fontsize=10, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 0].text(0.99, 0.85, label_count, transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            ax[idx, 0].text(0.99, 0.75, label_tail, transform=ax[idx, 0].transAxes, fontsize=9, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

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

    #add gridlines to all axes
    for row in range(ax.shape[0]):
        for col in range(ax.shape[1]):
            if (row, col) in axes_with_data:
                #ax[row, col].grid(axis='x',color='gray', linestyle='--', alpha=0.5)
                ax[row, col].set_facecolor('none')  # Remove plot background
            else:
                ax[row, col].grid(False)  # Disable grid for empty axes
    # Set x-ticks for specific axes and ensure they are visible
    for (row, col) in [(3,1), (4,0), (5,2)]:
        ax[row, col].set_xticks(np.arange(0, 300, 25))
        ax[row, col].tick_params(axis='x', which='both', length = 3, labelbottom=True)
        ax[row, col].set_xlabel('Frequency (Hz)', fontsize=12)
    #plt.subplots_adjust(hspace=0, wspace=0.1)  # Small vertical and horizontal space between subplots
    plt.xlim(5, 300)
    plt.tight_layout(pad=3, w_pad=0.5, h_pad=0)

    fig.savefig('histogram.png', dpi=300, bbox_inches='tight')
    plt.show()

