import numpy as np
import matplotlib.pyplot as plt

file_path = '/home/irseppi/REPOSITORIES/'
infile = open(file_path + 'parkshwynodal/input/node_crossings_db_UTM.txt', 'r')

equip_counts = {}  
flight_nums = {}   

for line in infile:
	data = line.split(',')  
	equip = data[-2]  # Get the equipment type from the line
	if equip == np.nan or equip == 'nan':
		equip = 'Unknown'
	flight_num = data[1]
	if flight_num not in flight_nums:
		flight_nums[flight_num] = []
		equip_counts[equip] = equip_counts.get(equip, 0) + 1  # move outside loop to get count of crossings instead of counts of flights

#Print summary statistics
print('In the crossings database:')
print('Number of equipment types:', len(equip_counts))
print('Number of flights:', len(flight_nums))

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

#Read in color text file to get different flights to be diffrent colors on map
colors=[]
with open(file_path + 'parkshwynodal/input/fig_style/colors.txt','r') as c_in:
	for i, line in enumerate(c_in):
		if (i + 1) % 9 == 0:
			c = str(line[0:7])
			colors.append(c)

# Plot the two pie charts side by side
fig, axes = plt.subplots(2, 1, figsize=(8, 14))

# Plot the first pie chart with the 'Other' category as the last slice
sorted_labels = list(equip_counts.keys())
sorted_labels.remove('Other')
sorted_labels.append('Other')  # Ensure 'Other' is the last slice

sorted_sizes = [equip_counts[label] for label in sorted_labels]

# Modify colors to make the 'Other' slice magenta
modified_colors = colors[4:(len(equip_counts)+4)]
modified_colors[-1] = 'magenta'  # Set the color for 'Other' slice to magenta

# Label wedges of pie chart with equipment type and flight number counts
wedges, _ = axes[0].pie(
	sorted_sizes,
	colors=modified_colors,
	textprops={'fontsize': 10},
)
# Manually move overlapping labels
for i, wedge in enumerate(wedges):
	if sorted_labels[i] == 'B77W':
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[0].text(
			x * 1.2, y * 1.05,
			f"{sorted_labels[i]}: {sorted_sizes[i]}",
			ha='center', va='center_baseline', fontsize=10,
		)
	elif sorted_labels[i] == 'DH3T':
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[0].text(
			x * 1.2, y * 1.11,
			f"{sorted_labels[i]}: {sorted_sizes[i]}",
			ha='center', va='center_baseline', fontsize=10,
		)
	elif sorted_labels[i] == 'DH3T':
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[0].text(
			x * 1.15, y * 1.1,
			f"{sorted_labels[i]}: {sorted_sizes[i]}",
			ha='center', va='center_baseline', fontsize=10,
		)
	elif sorted_labels[i] == 'AT73':
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[0].text(
			x * 1.5, y * 1.05,
			f"{sorted_labels[i]}: {sorted_sizes[i]}",
			ha='center', va='center_baseline', fontsize=10,
		)
	elif sorted_labels[i] == 'Unknown':
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[0].text(
			x * 1.4, y * 1.05,
			f"{sorted_labels[i]}: {sorted_sizes[i]}",
			ha='center', va='center_baseline', fontsize=10,
		)
	else:
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[0].text(
			x * 1.2, y * 1.07,
			f"{sorted_labels[i]}: {sorted_sizes[i]}",
			ha='center', va='center_baseline', fontsize=10,
		)
axes[0].axis('equal')
# Label wedges of pie chart with equipment type and flight number counts
wedges, _ = axes[1].pie(
	less_than_50.values(),
	colors=colors[(len(equip_counts)+10):((len(equip_counts)+10+len(less_than_50)))][::-1],
	textprops={'fontsize': 10},
)
# Manually move overlapping labels
for i, wedge in enumerate(wedges):
	label = list(less_than_50.keys())[i]
	size = list(less_than_50.values())[i]
	if label == 'B18T':
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[1].text(
			x * 1.2, y * 1.07,
			f"{label}: {size}",
			ha='center', va='center_baseline', fontsize=10,
		)
	elif label == 'R44':
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[1].text(
			x * 1, y * 1.09,
			f"{label}: {size}",
			ha='center', va='center_baseline', fontsize=10,
		)
	else:
		angle = (wedge.theta2 + wedge.theta1) / 2
		x = np.cos(np.deg2rad(angle))
		y = np.sin(np.deg2rad(angle))
		axes[1].text(
			x * 1.15, y * 1.07,
			f"{label}: {size}",
			ha='center', va='center_baseline', fontsize=10,
		)
axes[1].axis('equal')
plt.subplots_adjust(hspace=0.05)
plt.savefig(file_path + 'parkshwynodal/pie_charts_equipment.png', dpi=2500, bbox_inches='tight')

