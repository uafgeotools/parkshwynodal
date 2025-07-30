import os
import matplotlib.pyplot as plt
import numpy as np
#flights =  ['B763'] #['B737', 'B738', 'B739', 'B77W', 'B772', 'B788', 'B789', 'B763', 'B744','B733','B732','B77L','B748','CRJ2', 'A332', 'A359', 'E75S']
#flights = ['B190','BE20','C208','DH8A','AT73','SW4','PC12','DH3T','C441','B18T','B350','BE10']
#flights = ['AS50','R44'] 
#flights = ['PA31','DHC2','GA8','C180','C182','C206','C172','PA32','PA46','CH7B','PA30','C46','BE35','PA18','PA34'] #'C185'
flights = ['B737']

equip_overtone_dict = {}
count_dict = {}
color_dict = {}

y = 0
for eq in  flights:
    # Define the directory where your files are located
    file = eq + 'data_atmosphere_full.csv' #'output/Inversion_Results/'+eq+'data_atmosphere_full.txt'
    if eq not in equip_overtone_dict:
        equip_overtone_dict[eq] = []
        count_dict[eq] = []
        color_dict[eq] = np.random.rand(3,)
    with open(file, 'r') as f:

        # Read the data from the file and append it to the list
        for line in f.readlines():
            y += 1
            data = []
            counts = []
            lines = line.split(',')
            peaks = np.array(lines[9])
            peaks = str(peaks)  # Replace "string" with "str"
            peaks = np.array(peaks.split(' '))
            for peak in peaks:
                try:
                    peak = float(peak)
                except:
                    continue
                data.append(peak)
                counts.append(float(y))
            equip_overtone_dict[eq].extend(data)
            count_dict[eq].extend(counts)

plt.figure(figsize=(10, 6))
for equip, peaks in equip_overtone_dict.items():
    x = np.array(peaks)
    y = np.array(count_dict[equip])
    color = color_dict[equip]
    plt.scatter(peaks, y, c=color)
plt.legend(equip_overtone_dict.keys(), loc='upper right', fontsize='small')
plt.xlim(0, 305)
plt.xticks(np.arange(0, 305, 5))
plt.grid(True)
plt.show()


plt.figure()
for equip, peaks in equip_overtone_dict.items():
    x = np.array(peaks)

    color = color_dict[equip]
    plt.hist(peaks, color=color, bins=100)
plt.legend(equip_overtone_dict.keys(), loc='upper right', fontsize='small')
plt.xticks(np.arange(0, 305, 5))
plt.xlim(0,305)
#plt.ylim(0.7, 7)
plt.show()
