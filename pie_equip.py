import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

infile = open('input/node_crossings_db_UTM.txt', 'r')
#outfile = open('output.csv', 'a')  # Open the file in append mode

equip_counts = {}  # Define the "equip_counts" dictionary before the loop
flight_nums = {}  # Define the "equip_counts" dictionary before the loop

for line in infile:
	data = line.split(',')  # Split the line by commas
	equip = data[-2]  # Get the equipment type from the line
	if equip == np.nan or equip == 'nan':
		equip = 'Unknown'
	flight_num = data[1]
	if flight_num not in flight_nums:
		flight_nums[flight_num] = []
		#flight_nums[equip].extend([flight_num]) 
		equip_counts[equip] = equip_counts.get(equip, 0) + 1  # move outside loop to get count of crossings instead of counts of flights
#equip_counts = {'Unkown': 296, 'AT73': 50, 'B190': 130, 'B738': 136, 'B737': 367, 'B739': 139, 'C185': 61, 'B77W': 67, 'PA31': 176, 'DH8A': 95, 'DHC6': 1, 'C172': 13, 'C208': 244, 'DH3T': 40, 'BE20': 37, 'E75S': 4, 'AS50': 4, 'B744': 9, 'SW4': 27, 'PC12': 52, 'B772': 33, 'B789': 30, 'DHC2': 12, 'B407': 12, 'R44': 5, 'A359': 10, 'B06': 2, 'C46': 2, 'B763': 9, 'GA8': 20, 'B732': 1, 'C182': 13, 'CH7B': 5, 'B788': 18, 'C441': 7, 'B18T': 4, 'C180': 8, 'PA34': 1, 'B77L': 6, 'B350': 1, 'PA18': 22, 'C206': 7, 'BE35': 1, 'C82S': 1, 'B733': 6, 'PA46': 3, 'A332': 1, 'PA32': 9, 'GLF5': 1, 'B748': 1, 'CRJ2': 1, 'AT8T': 3, 'BE10': 1, 'AC6L': 4, 'B412': 2, 'PA30': 2, 'BE58': 1, 'BE36': 1}
print(equip_counts)
print(len(equip_counts))
print(len(flight_nums))
# Plotting the first pie chart
labels = equip_counts.keys()
sizes = equip_counts.values()

# Create a new dictionary for values less than 50
less_than_50 = {label: size for label, size in zip(labels, sizes) if size < 10}
less_than_50 = {t: g for t, g in sorted(less_than_50.items(), key=lambda item: item[1], reverse=True)}

# Calculate the sum of values less than 50
less_than_50_sum = sum(less_than_50.values())

# Add the 'Other' category to the dictionary with the sum of values less than 50
equip_counts['Other'] = less_than_50_sum

# Remove the keys with values less than 15 from the original dictionary
equip_counts = {label: size for label, size in equip_counts.items() if size >= 10}
equip_counts = {k: v for k, v in sorted(equip_counts.items(), key=lambda item: item[1], reverse=True)}
# Define a color dictionary
colors=[]
#Read in color text file to get different flights to be diffrent colors on map
with open('input/colors.txt','r') as c_in:
	for i, line in enumerate(c_in):
		if (i + 1) % 9 == 0:
			c = str(line[0:7])
			colors.append(c)

# Plot the two pie charts side by side
fig, axes = plt.subplots(1, 2, figsize=(24, 12))

# Plot the first pie chart with the 'Other' category as the last slice
sorted_labels = list(equip_counts.keys())
sorted_labels.remove('Other')
sorted_labels.append('Other')  # Ensure 'Other' is the last slice

sorted_sizes = [equip_counts[label] for label in sorted_labels]

axes[0].pie(sorted_sizes, labels=[f"{label}: {size}" for label, size in zip(sorted_labels, sorted_sizes)], colors=colors[4:(len(equip_counts)+4)])
axes[0].set_title("Equipment Counts (Including 'Other')")

# Plot the second pie chart for values less than 50
axes[1].pie(less_than_50.values(), labels=[f"{label}: {size}" for label, size in less_than_50.items()], colors=colors[(len(equip_counts)+10):((len(equip_counts)+10+len(less_than_50)))][::-1])
axes[1].set_title("Equipment Counts (Less than 50)")

plt.tight_layout()
plt.show()
