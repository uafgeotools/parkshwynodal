import pandas as pd
import numpy as np

infile = open('input/node_crossings_db_UTM.txt', 'r')

equip_counts = {}  # Define the "equip_counts" dictionary before the loop
flight_nums = {}  # Define the "flight_nums" dictionary before the loop
crossings = {}
engine_types = {}
engine_counts = {}
wake_turbs = {}

for line in infile:
	data = line.split(',')  # Split the line by commas
	equip = data[-2]  # Get the equipment type from the line
	if equip not in equip_counts:
		file = open('output/inv_results_old/' + equip + '_full_inv_results.csv', 'r')
		for ll in file:
			l = ll.split(',')
			flight_num = l[0]
			sta = l[2]
			closest_time = l[3]
			key = (flight_num, sta, closest_time)
			crossings[key] = crossings.get(key, 0) + 1
	if equip == np.nan or equip == 'nan':
		equip = 'Unknown'
	flight_num = data[1]
	if flight_num not in flight_nums:
		flight_nums[equip] = []
		flight_nums[equip].extend([flight_num]) 
	equip_counts[equip] = equip_counts.get(equip, 0) + 1  # move outside loop to get count of crossings instead of counts of flights

print(equip_counts)
print(len(equip_counts))
print(len(flight_nums))



