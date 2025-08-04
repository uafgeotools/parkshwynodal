import matplotlib.pyplot as plt
import numpy as np
flights = ['PA31','DHC2','GA8','C180','C182','C206','C172','PA32','PA46','CH7B','PA30','C185','AS50','R44','B190','BE20','C208','DH8A','AT73','SW4','PC12','DH3T','C441','B18T','BE10','B763','B737', 'B738', 'B739', 'B763','B733']

equip_overtone_dict = {}

for eq in flights:
    # Define the directory where your files are located
    file = 'output/inv_results/' + eq + '_full_inv_results.csv' 
    if eq not in equip_overtone_dict:
        equip_overtone_dict[eq] = []
    with open(file, 'r') as f:
        # Read the data from the file and append it to the list
        data = []
        for line in f.readlines():
            lines = line.split(',')

            if lines[-2] == "Forward Model":
                continue
            peaks = np.array(lines[9])
            peaks = str(peaks)  # Replace "string" with "str"
            # remove[ and ] from the string
            peaks = peaks.replace('[', '').replace(']', '')
            peaks = np.array(peaks.split(' '))
            for peak in peaks:
                if peak == '':
                    continue
                peak = float(peak)
            
                data.append(peak)
        equip_overtone_dict[eq].extend(data)


print(len(equip_overtone_dict.keys()), "equipments found")
fig, ax = plt.subplots(int(len(equip_overtone_dict.keys())/3), 3, figsize=(20, 25), sharex=True)

for i, (equip, peaks) in enumerate(equip_overtone_dict.items()):
    bins = np.arange(min(peaks), max(peaks) + 3, 3)
    ax[i//3, i%3].hist(peaks, color='k', bins=bins, alpha=0.5, edgecolor='black')
    ax[i//3, i%3].text(0.01, 0.95, equip, transform=ax[i//3, i%3].transAxes, fontsize=10, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    counts, _ = np.histogram(peaks, bins=bins)
    ax[i//3, i%3].set_yticks([counts.max()])
plt.subplots_adjust(hspace=0, wspace=0.1)  # Small vertical and horizontal space between subplots
plt.xlim(5,300)
fig.savefig('histogram.png', dpi=300, bbox_inches='tight')
plt.show()

