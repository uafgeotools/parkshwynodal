import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patheffects as patheffects

file = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/output/NGT_flight_param_inv_DB.txt', sep=",")
flight_nums = file['Flight_Number']
freq_peaks = file['Meas_Source_Frequency_Array']
varr_dict = file['Variance']

flight_num_dict = {'10512184': ['527958214', '529754214', '529970458', '530342801', '530489496', '530681886', '530893864', 
                                '531043310', '531236830', '531254054', '533770152', '533950817', 
                                '534899218'], 
                    '10572742': ['528502194', '528518474', '528510459', '528711629', '528698927', '529409728', 
                                '529416700', '530144820', '530501643', '530672521', '530695031', '530706111', 
                                '530842926', '531272879', '531417015', '531424585', '531583949', 
                                '531605202', '531901801', '531774192', '533421550']}

# Create a dictionary to store the color for each tail number
color_dict = {}
peaks_dict = {}
all_med = {}
error_dict = {}
sort_1 = []
sort_2 = []
sort_3 = []

# Iterate over each line in the file
for i in range(len(flight_nums)):
    flight_num = str(flight_nums[i])
    if flight_num in flight_num_dict['10512184']:
        tail_num = '10512184'
    elif flight_num in flight_num_dict['10572742']:
        tail_num = '10572742'
    else:
        continue
    if tail_num not in color_dict:
        color_dict[tail_num] = []
        peaks_dict[tail_num] = []
        all_med[tail_num] = []
        error_dict[tail_num] = []

    peaks = np.array(freq_peaks[i].strip('[]').split(' '), dtype=float)
    Cpost0 = np.array(varr_dict[i].strip('[]').split(' '), dtype=float)

    f1 = []
    peak_old = 0
    for tt, peak in enumerate(peaks):
        if len(peaks) == 0 or peak == peaks[0]:
            peak_old = float(peak)
            continue
        diff = float(peak) - float(peak_old)
        f1.append(diff)
        peak_old = float(peak)
    peaks_dict[tail_num].extend(peaks.tolist())
    all_med[tail_num].extend([np.nanmedian(f1)])
    #Generate random samples of f0 values withing their sigma from the covariance matrix 
    #Calculate the median of the differences and MAD to obtain error
    f_range = []
    NTRY = 1000
    for N in range(NTRY):
        ftry = []
        for c_index  in range(4, len(Cpost0)):
            xmin = peaks[c_index-4] - Cpost0[c_index]
            xmax = peaks[c_index-4] + Cpost0[c_index]
            xtry = xmin + (xmax-xmin)*np.random.rand()
            ftry.append(xtry)

        ftry = np.sort(ftry)
        f1 = []
        for g in range(len(ftry)):
            if g == 0:
                continue
            diff = ftry[g] - ftry[g - 1]
            f1.append(diff)
        f_range.append(np.nanmedian(f1))
    med_df = np.nanmedian(f_range)
    mad_df = np.nanmedian(np.abs(f_range - med_df))

    error_dict[tail_num].extend([mad_df])

    if str(tail_num) == '10512184' and med_df < 20:
        sort_1.append(mad_df)
    elif str(tail_num) == '10512184' and med_df > 20:
        sort_2.append(mad_df)
    elif str(tail_num) == '10572742':
        sort_3.append(mad_df)

print('Count for 10572742:', len(peaks_dict['10572742']))
print('Count for 10512184:', len(peaks_dict['10512184']))
print('Flight num for 10572742:', len(flight_num_dict['10572742']))
print('Flight num for 10512184:', len(flight_num_dict['10512184']))
print('Crossings for 10512184:', len(sort_1+sort_2))
print('Crossings for 10572742:', len(sort_3))

fig,ax1 = plt.subplots(1, 1, sharex=False, figsize = (50,20)) #figsize=(50,20))     

ax1.margins(x=0)
ax2 = fig.add_axes([0.87, 0.072, 0.125, 0.904], sharey=ax1)

color_dict['10512184'] = [1.0, 0.5, 0.0]  # Orange color in RGB
color_dict['10572742'] = [0.0, 0.5, 1.0]  # Blue color in RGB

for tail_num, peaks in peaks_dict.items():
    error_med = np.nanmedian(error_dict[tail_num])
    if str(tail_num) == '10512184':
        sort_1 = np.nanmedian(np.array(sort_1))
        sort_2 = np.nanmedian(np.array(sort_2))
        print(f'Tail Number: {tail_num}, Median Error: {error_med}, Sort 1: {sort_1}, Sort 2: {sort_2}')
    elif str(tail_num) == '10572742':
        print(f'Tail Number: {tail_num}, Median Error: {error_med}')
    color = color_dict[tail_num]
    med = all_med[tail_num]

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
