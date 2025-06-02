import os
import matplotlib.pyplot as plt
import numpy as np


# Define the directory where your files are located
file = 'C185data_atmosphere.csv'

# Initialize a list to store the data from the files
data = []
with open(file, 'r') as f:
    x = 0
    # Read the data from the file and append it to the list
    for line in f.readlines():
        lines = line.split(',')
        #quality_num = int(lines[8]) #lines[9] for output2

        #if quality_num < 21:
        #    continue
        try:

            peaks = np.array(lines[7])
            
            
        except:
            continue
        peaks = str(peaks)  # Replace "string" with "str"
        peaks = np.array(peaks.split(' '))
        for peak in peaks:
            try:
                peak = float(peak)
            except:
                continue
            print(peak)
            data.append(peak)
        x += 1
plt.figure()
# Plot the histogram
plt.hist(np.round(data,1),250)

# Add tick marks at every value of 5
plt.xticks(np.arange(15, 250, 5))

plt.show()
