import os
import matplotlib.pyplot as plt
import numpy as np

# Define the directory where your files are located
file = 'C185data_atmosphere.csv'
plt.figure()
# Initialize lists to store the data from the files
mode1_x = []
mode1_y = []
mode2_x = []
mode2_y = []
mode3_x = []
mode3_y = []
mode1_f1 = []
mode1_f2 = []
mode1_f3 = []
mode1_f4 = []
mode1_f5 = []
mode1_f6 = []
mode1_f7 = []
mode1_f8 = []
mode1_f9 = []
mode1_f10 = []
mode1_f11 = []
mode2_f1 = []
mode2_f2 = []
mode2_f3 = []
mode2_f4 = []
mode2_f5 = []
mode2_f6 = []
mode2_f7 = []
mode2_f8 = []
mode2_f9 = []
mode2_f10 = []
mode2_f11 = []
mode2_f12 = []
mode2_f13 = []
mode2_f14 = []
QNUM = False
with open(file, 'r') as f:
    x = 0
    xx = 0
    yy = 0
    zz = 0
    # Read the data from the file and append it to the respective lists based on the mode
    for line in f.readlines():
        lines = line.split(',')
        if QNUM == True:
            quality_num = int(lines[8]) #lines[9] for output2
            if int(quality_num) < 21:
                continue
            else:
                print(quality_num)
        x += 1
        mode_x = 0
        peaks = np.array(lines[7])

        peaks = str(peaks)  

        peaks = np.char.replace(peaks, '[', '')
        peaks = np.char.replace(peaks, ']', '')
        peaks = str(peaks)

        peaks = np.array(peaks.split(' '))
        for peak in peaks:
            if peak == '' or peak == ' ' or peak == '   ':
                continue
            try:
                peak = float(peak[0:-1])
            except:
                continue
            if peak < 111 or peak > 126:
                continue
            else:
                if abs(float(peak) - 124.5) <= abs(float(peak) - 114):
                    xx += 1
                    mode_x = 1
                    break
                elif abs(float(peak) - 114) <= abs(float(peak) - 124.5):
                    mode_x = 2
                    yy += 1
                    break
                else:
                    mode_x = 3
                    continue
                    zz += 1
            
        if mode_x == 1:
            for peak in peaks:
                try:
                    peak = float(peak[0:-1])
                except:
                    continue
                if abs(float(peak) - 62) <= 1.5:
                    mode1_f1.append(peak)
                elif abs(float(peak) - 82.5) <= 1.5:
                    mode1_f2.append(peak)
                elif abs(float(peak) - 103) <= 2:
                    mode1_f3.append(peak)
                elif abs(float(peak) - 124.5) <= 1.5:
                    mode1_f4.append(peak)
                elif abs(float(peak) - 145) <= 2:
                    mode1_f5.append(peak)
                elif abs(float(peak) - 166) <= 2:
                    mode1_f6.append(peak)
                elif abs(float(peak) - 186.5) <= 2:
                    mode1_f7.append(peak)
                elif abs(float(peak) - 208) <= 2:
                    mode1_f8.append(peak)
                elif abs(float(peak) - 228.5) <= 2:
                    mode1_f9.append(peak)
                elif abs(float(peak) - 249) <= 4:
                    mode1_f10.append(peak)
                elif abs(float(peak) - 269) <= 4:
                    mode1_f11.append(peak)
                mode1_x.append(peak)
                mode1_y.append(xx)
        elif mode_x == 2:
            for peak in peaks:
                try:
                    peak = float(peak[0:-1])
                   
                except:
                    continue
                mode2_x.append(peak)
                mode2_y.append(yy)
                
                if abs(float(peak) - 23) <= 4:
                    mode2_f1.append(peak)
                elif abs(float(peak) - 39) <= 4:
                    mode2_f2.append(peak)
                elif abs(float(peak) - 58) <= 4:
                    mode2_f3.append(peak)
                elif abs(float(peak) - 77) <= 4:
                    mode2_f4.append(peak)
                elif abs(float(peak) - 96.5) <= 4:
                    mode2_f5.append(peak)
                elif abs(float(peak) - 116) <= 5:
                    mode2_f6.append(peak)
                elif abs(float(peak) - 133) <= 6:
                    mode2_f7.append(peak)
                elif abs(float(peak) - 153.5) <= 4:
                    mode2_f8.append(peak)
                elif abs(float(peak) - 174) <= 3:
                    mode2_f9.append(peak)
                elif abs(float(peak) - 192) <= 6:
                    mode2_f10.append(peak)
                elif abs(float(peak) - 215) <= 6:
                    mode2_f11.append(peak)
                elif abs(float(peak) - 233) <= 6:
                    mode2_f12.append(peak)
                elif abs(float(peak) - 255) <= 6:
                    mode2_f13.append(peak)
        else:
            for peak in peaks:
                try:
                    peak = float(peak)
                except:
                    continue
                mode3_x.append(peak)
                mode3_y.append(zz)


medians_mode1 = [np.nanmedian(mode1_f1), np.nanmedian(mode1_f2), np.nanmedian(mode1_f3), np.nanmedian(mode1_f4), np.nanmedian(mode1_f5), np.nanmedian(mode1_f6), np.nanmedian(mode1_f7), np.nanmedian(mode1_f8), np.nanmedian(mode1_f9), np.nanmedian(mode1_f10), np.nanmedian(mode1_f11)]
medians_mode2 = [np.nanmedian(mode2_f1), np.nanmedian(mode2_f2), np.nanmedian(mode2_f3), np.nanmedian(mode2_f4), np.nanmedian(mode2_f5), np.nanmedian(mode2_f6), np.nanmedian(mode2_f7), np.nanmedian(mode2_f8), np.nanmedian(mode2_f9), np.nanmedian(mode2_f10), np.nanmedian(mode2_f11), np.nanmedian(mode2_f12), np.nanmedian(mode2_f13), np.nanmedian(mode1_f11)]
medians_mode1_rounded = [round(median, 1) for median in medians_mode1 if not np.isnan(median)]

medians_mode2_rounded = [round(median, 1) for median in medians_mode2 if not np.isnan(median)]

for medians in medians_mode1:
    plt.axvline(x=medians, color='orange', linestyle='--')
for medians in medians_mode2:
    plt.axvline(x=medians, color='blue', linestyle='--')

# Set the x-axis labels as rounded median values
plt.xticks(medians_mode1_rounded + medians_mode2_rounded)

# Add labels for the differences between medians
for i in range(len(medians_mode1_rounded) - 1):
    diff = round(medians_mode1_rounded[i+1] - medians_mode1_rounded[i], 1)
    plt.text(medians_mode1_rounded[i]+(diff/2), 1, f'{diff}', ha='center', va='bottom')
    plt.annotate('', xy=(medians_mode1_rounded[i], 1), xytext=(medians_mode1_rounded[i+1], 1),
                 arrowprops=dict(arrowstyle='<->', color='orange'))
    
for i in range(len(medians_mode2_rounded) - 1):
    diff = round(medians_mode2_rounded[i+1] - medians_mode2_rounded[i], 1)
    plt.text(medians_mode2_rounded[i]+(diff/2), 3, f'{diff}', ha='center', va='bottom')
    plt.annotate('', xy=(medians_mode2_rounded[i], 3), xytext=(medians_mode2_rounded[i+1], 3),
                 arrowprops=dict(arrowstyle='<->', color='blue'))
diff_mode1_mode2 = False
if diff_mode1_mode2 == True:
    for i in range(0,7):
        diff = round(medians_mode1_rounded[i] - medians_mode2_rounded[i+2], 1)
        plt.text(medians_mode1_rounded[i]-(diff/2), 5, f'{diff}', ha='center', va='bottom')
        plt.annotate('', xy=(medians_mode1_rounded[i], 5), xytext=(medians_mode2_rounded[i+2], 5),
                    arrowprops=dict(arrowstyle='<->', color='green'))

# Plot the dots for modes
plt.scatter(mode1_x, mode1_y, label='Mode 1', color='orange')
mode2_y_modified = [(np.max(mode1_y)) + y for y in mode2_y]
plt.scatter(mode2_x, mode2_y_modified, label='Mode 2', color='blue')

# Add legend
plt.legend()


plt.show()
