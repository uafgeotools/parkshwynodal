import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
file = open('old_new_invers.txt','r')

color_dict = {}
count_dict = {}
peaks_dict_new = {}
y = 0
#plt.figure()
# Iterate over each line in the file
for line in file.readlines():
    y += 1
    counts = []

    pppp_new = []

    lines = line.split(',')
    flight_num = float(lines[1])
    peaks_new = np.array(lines[12])

    print(peaks_new)
    y += 1
    counts = []
    pppp_new = []
    if flight_num not in color_dict:
        color_dict[flight_num] = np.random.rand(3,)
        peaks_dict_new[flight_num] = []
        count_dict[flight_num] = []

    peaks = str(peaks_new)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))
    
    for j in range(len(peaks)):
        peak_new = peaks[j]
        pppp_new.append(float(peak_new))
        counts.append(float(y))
    peaks_dict_new[flight_num].extend(pppp_new)  # Convert peaks_new[i] to a list
    count_dict[flight_num].extend(counts)  # Create a list of the same length as peaks_new[i]

fig, ax1 = plt.subplots(1, 1, sharex=False, figsize=(8, 6))
#ax1.margins(x=0)
ax2 = fig.add_axes([0.83, 0.11, 0.07, 0.77], sharey=ax1)
ax3 = fig.add_axes([0.90, 0.11, 0.05, 0.77], sharey=ax1)
ax4 = fig.add_axes([0.95, 0.11, 0.04, 0.77], sharey=ax1)
ax1.set_title('Frequency Peaks')
difff = []
for flight_num, peaks in peaks_dict_new.items():
    color = color_dict[flight_num]
    y = count_dict[flight_num]

    peaks = np.array(peaks)  # Convert peaks to a NumPy array
    y = np.array(y)  # Convert y to a NumPy array

    ax1.scatter(peaks, y, c=color, label=flight_num)
    f1 = []
    for i,peak in enumerate(peaks):
        if peak < 111 or peak > 126:
                continue
        else:
            rpm = 60 * (peak/3)
            ax2.scatter(rpm, y[i], c=color)
        if y[i] == y[i-1]:
            diff = float(peak) - float(peaks[i-1])

            f1.append(diff)

            if not np.isnan(np.nanmedian(f1)):
                del_f = (np.nanmedian(f1))
                if del_f > 24 or del_f < 16:
                    continue
                ax3.scatter(del_f, y[i], c=color)
                ax4.scatter(peak/del_f, y[i], c=color)
                difff.append(peak/del_f)
        else:
            f1 = []
            continue

ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax2.set_title('rpm')
ax3.set_title('\u0394' + 'F')
ax1.set_xlabel('Frequency')
ax2.set_xlabel('rpm')
ax3.set_xlabel('\u0394' + 'F')
ax4.set_xlabel('F0 /'+ '\u0394' + 'F')
ax1.legend(loc='upper left', fontsize='x-small')
ax1.set_xlim(0, 300)
ax1.set_xticks(range(0, 251, 25))

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
plt.show()
print(len(difff), len(counts))
y = counts
plt.figure()
for i in range(len(difff)):
    plt.scatter(np.median(np.absolute(difff[i] - np.median(difff))), y[i], c='b')
plt.show()