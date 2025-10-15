import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patheffects as patheffects

file = open('/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results_ngt/C185_full_inv_results.txt', 'r')

file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

# Create a dictionary to store the color for each tail number
color_dict = {}
peaks_dict = {}
all_med = {}
flight_num_hold = {}
error_dict = {}
sort_1 = []
sort_2 = []
count_id1 = 0
count_id2 = 0
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]
    if lines[11] == '00':
        continue
    peaks = np.array(lines[9])
    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')
    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    error_strs = [e for e in lines[10].strip('[]').split(' ') if e.strip() != '']
    Cpost0 = np.array(error_strs[3:], dtype=float)

    ppp = []
    f1 = []
    peak_old = 0
    for tt, peak in enumerate(peaks):
        if np.abs(float(peak) - float(peak_old))< 10:
            continue
        ppp.append(float(peak))

        if len(peaks) == 0 or peak == peaks[0]:
            peak_old = float(peak)
            continue

        diff = float(peak) - float(peak_old)
        f1.append(diff)
        peak_old = float(peak)
    #Generate random samples of f0 values withing their sigma from the covariance matrix 
    #Calculate the median of the differences and MAD to obtain error
    f_range = []

    NTRY = 1000
    for N in range(NTRY):
        ftry = []
        for c_index  in range(4, len(Cpost0)):
            xmin = ppp[c_index-4] - (Cpost0[c_index]/len(error_strs))
            xmax = ppp[c_index-4] + (Cpost0[c_index]/len(error_strs))
            xtry = xmin + (xmax-xmin)*np.random.rand()
            ftry.append(xtry)

        ftry = np.sort(ftry)
        f_hold = []
        for g in range(len(ftry)):
            if g == 0:
                continue
            diff = ftry[g] - ftry[g - 1]
            f_hold.append(diff)
        med = np.nanmedian(f_hold)
        f_range.append(med)
    med_df = np.nanmedian(f_range)
    mad_df = np.nanmedian(np.abs(f_range - med_df))

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = []
                peaks_dict[tail_num] = []
                all_med[tail_num] = []
                flight_num_hold[tail_num] = []
                error_dict[tail_num] = []
                break
        else:
            continue
    if tail_num == '10572742':
        count_id1 += 1
    elif tail_num == '10512184':
        count_id2 += 1
        
    peaks_dict[tail_num].extend(ppp)
    all_med[tail_num].extend([np.nanmedian(f1)])
    error_dict[tail_num].extend([mad_df])
    if flight_num not in flight_num_hold[tail_num]:
        flight_num_hold[tail_num].append(flight_num)
    if str(tail_num) == '10512184' and med_df < 20:
        sort_1.append(mad_df)
    elif str(tail_num) == '10512184' and med_df > 20:
        sort_2.append(mad_df)
print('Count for 10572742:', count_id1)
print('Count for 10512184:', count_id2)
fig,ax1 = plt.subplots(1, 1, sharex=False, figsize = (50,20)) #figsize=(50,20))     

ax1.margins(x=0)
ax2 = fig.add_axes([0.87, 0.072, 0.125, 0.904], sharey=ax1)

pos = 1
tail_num_hold = 0
color_dict[10512184] = [1.0, 0.5, 0.0]  # Orange color in RGB 
color_dict[10572742] = [0.0, 0.5, 1.0]  # Blue color in RGB

for tail_num, peaks in peaks_dict.items():
    error_med = np.nanmedian(error_dict[tail_num])
    if str(tail_num) == '10512184':
        print(sort_1, sort_2)
        sort_1 = np.nanmedian(np.array(sort_1))
        sort_2 = np.nanmedian(np.array(sort_2))
        print(f'Tail Number: {tail_num}, Median Error: {error_med}, Sort 1: {sort_1}, Sort 2: {sort_2}')
    print(f'Tail Number: {tail_num}, Median Error: {error_med}')
    color = color_dict[tail_num]
    med = all_med[tail_num]
    if str(tail_num) != '10572742' and str(tail_num) != '10512184':
        continue
    ax1.hist(peaks, bins=270, color=color, alpha=0.8, label=tail_num, zorder = 10)  
    ax2.hist(med, bins=270, color=color, alpha=0.8, zorder = 10)  
    ax1.hist(peaks, bins=270, color=color, histtype='step',zorder = 15)  
    ax2.hist(med, bins=270, color=color, histtype='step', zorder = 15)  

ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax1.set_xlabel('Frequency (Hz)', fontsize=30)
ax2.set_xlabel('Median '+'\u0394'+'f (Hz)', fontsize=30)
ax1.set_xticks(range(10, 270, 10))
ax1.set_xticklabels(
    [str(x) if x % 20 == 0 else '' for x in range(10, 270, 10)],
    fontsize=30
)
ax1.set_yticks(range(0, 90, 20))
ax1.tick_params(axis='both', labelsize=30)  # Increase font size for tick labels
ax2.tick_params(axis='both', labelsize=30)  # Increase font size for tick labels
ax2.set_xticks(np.arange(18.5, 22, 0.5))
ax2.set_xticklabels(
    [str(int(x)) if x % 1 == 0 else '' for x in np.arange(18.5, 22, 0.5)],
    fontsize=30
)
ax2.set_xlim(18, 22)
ax1.set_ylim(0, 82)

del_f_t1 = 19.6
del_f_t2_1 = 19.1 #2
del_f_t2_2 = 20.55 #57
x_label = []
for g in range(0,14):
    ax1.axvline(x= (1 + g) * del_f_t1, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
    if g != 13:
        x_label.append((1 + g) * del_f_t1)
    ax1.axvline(x= (1 + g) * del_f_t2_1, color = [1.0, 0.5, 0.0], ls = '--', zorder=0, linewidth=1)
    ax1.axvline(x= (1 + g) * del_f_t2_2, color = [1.0, 0.5, 0.0], ls = '--', zorder=0, linewidth=1)
# Add labels inside the plot at the top for each overtone position
for x in x_label:
    ax1.text(
        x,
        ax1.get_ylim()[1] - 1,
        f"{x:.1f}",
        color=[0.0, 0.5, 1.0],
        fontsize=30,
        ha='center',
        va='top',
        rotation=0,
        path_effects=[patheffects.withStroke(linewidth=2, foreground='white')],
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2')
    )

for label in ax1.get_yticklabels():
    label.set_fontsize(30)
ax1.set_xlim(5,310)
med_label = 19.6
ax2.text(
    med_label,
    ax2.get_ylim()[1] - 1,  # Place at the very top
    f"{med_label:.1f}",
    color=[0.0, 0.5, 1.0],
    fontsize=30,
    ha='center',
    va='top',
    zorder=100,  # Ensure it's on top of other plot elements
    rotation=0,
    path_effects=[patheffects.withStroke(linewidth=2, foreground='white')],
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2')
)
# Change tick color, outline in black, make them bold, and increase font size
for label in ax1.get_xticklabels():
    label.set_fontsize(30)

for label in ax2.get_xticklabels():
    label.set_fontsize(30)
print(del_f_t1, del_f_t2_1, del_f_t2_2)
ax2.axvline(x=del_f_t1, color = [0.0, 0.5, 1.0], ls = '--', zorder=0, linewidth=1)
ax2.axvline(x=del_f_t2_1, color =   [1.0, 0.5, 0.0], ls = '--', zorder=0, linewidth=1)
ax2.axvline(x=del_f_t2_2, color =  [1.0, 0.5, 0.0], ls = '--', zorder=0, linewidth=1)
plt.tight_layout(pad=3.5, w_pad=0.5, h_pad=1.5)
plt.show()
