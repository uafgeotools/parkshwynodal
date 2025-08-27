import pandas as pd
import os
jets = True
with open('input/node_crossings_db_UTM.txt', 'r') as infile:
	infile_lines = infile.readlines()

Equipment = ['DH8A', 'B190','BE20','PC12','DH3T','C208','AT73','SW4','C441','B18T','B350','BE10','AS50','R44','C185','PA31','DHC2','GA8','C180','C182','C206','C172','PA32','PA46','CH7B','PA30','C46','BE35','PA18','PA34']
Equipment_jets = ['B737','B738','B739','B77W','B772','B789','B788','B733','B763','A359','B77L','B744','E75S','B732','A332','B748']

equip_data = pd.read_csv('input/20231010_Aircraft_UA_Fairbanks.csv', sep=",")
man = equip_data['MANUFACTURER']
model = equip_data['Model']
engine = equip_data['Engine Type']
engine_count = equip_data['Engine Count']
wake_turb = equip_data['Wake Turbulence Category']
equip_type = equip_data['Type Designator']

engine_counts = {}
wake_turbs = {}
type_name = {}
engine_type_dict = {}
man_dict = {}

flight_nums = {}  # Define the "flight_nums" dictionary before the loop
crossings = {}
tail_numbers = {}
nodes_count = {}

flight_nums_inverted = {}
crossings_inverted = {}
tail_numbers_inverted = {}
node_inverted = {}
if not jets:
	for key in Equipment:
		flight_nums[key] = []
		crossings[key] = 0
		tail_numbers[key] = []
		nodes_count[key] = []
		flight_nums_inverted[key] = []
		crossings_inverted[key] = 0
		tail_numbers_inverted[key] = []
		node_inverted[key] = []
		for eq_type in equip_type:
			if str(eq_type) == str(key):
				ind = equip_type[equip_type == eq_type].index[0]
				engine_counts[key] = engine_count[ind]
				wake_turbs[key] = wake_turb[ind]
				type_name[key] = model[ind]
				engine_type_dict[key] = engine[ind]
				man_dict[key] = man[ind]
				break
			else:
				continue
		for line in infile_lines:
			data = line.split(',')  # Split the line by commas
			equip = data[10]  # Get the equipment type from the line
			if equip not in Equipment:
				continue
			if equip == key:
				crossings[key] = crossings.get(key, 0) + 1
				flight_num = data[1]
				if flight_num not in flight_nums[key]:
					flight_nums[key].extend([flight_num])
				sta = data[9]
				date = data[0]
				flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_flights.csv', sep=",")
				flight = flight_data['flight_id']
				tailnumber = flight_data['aircraft_id']
				for i, f_id in enumerate(flight):
					if str(f_id) == str(flight_num):
						tailnum = tailnumber[i]
						break
					else:
						continue
				if tailnum not in tail_numbers[key]:
					tail_numbers[key].extend([tailnum])
				if sta not in nodes_count[key]:
					nodes_count[key].extend([sta])

		inv_file = 'output/inv_results/' + key + '_full_inv_results.txt'
		if not os.path.exists(inv_file):
			continue
		file = open('output/inv_results/' + key + '_full_inv_results.txt', 'r')
		for ll in file:
			l = ll.split(',')
			if str(l[-1]) == 'Forward Model':
				continue
			crossings_inverted[key] = crossings_inverted.get(key, 0) + 1
			f_num = l[1]
			if f_num not in flight_nums_inverted[key]:
				flight_nums_inverted[key].extend([f_num])
			sta = l[2]
			if sta not in node_inverted[key]:
				node_inverted[key].extend([sta])
			date = l[0]
			flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_flights.csv', sep=",")
			flight = flight_data['flight_id']
			tailnumber = flight_data['aircraft_id']
			for i, f_id in enumerate(flight):
				if str(f_id) == str(f_num):
					tailnum = tailnumber[i]
					break
				else:
					continue
			if tailnum not in tail_numbers_inverted[key]:
				tail_numbers_inverted[key].extend([tailnum])

	final_table = open('paper_table.txt', 'w')
	# Write the header row
	final_table.write('Equipment,Manufacturer,Model,Engine Type,Engine Count,Wake Turbulence Category,Total Flights,Inverted Flights,Total Crossings,Inverted Crossings,Total Tail Numbers,Inverted Tail Numbers,Total Nodes,Inverted Nodes\n')
	for eq in Equipment:

		final_table.write(eq +',' + str(man_dict[eq]) + ',' + str(type_name[eq]) + ',' + str(engine_type_dict[eq]) + ',' + str(engine_counts[eq]) + ',' + str(wake_turbs[eq]) + ',' + str(len(flight_nums[eq])) + ',' + str(len(flight_nums_inverted[eq])) + ',' + str(crossings[eq]) + ',' + str(crossings_inverted[eq]) + ',' + str(len(tail_numbers[eq])) + ',' + str(len(tail_numbers_inverted[eq])) + ',' + str(len(nodes_count[eq])) +  ',' + str(len(node_inverted[eq])) + '\n')

	final_table.close()
if jets:
	for key in Equipment_jets:
		flight_nums[key] = []
		crossings[key] = 0
		tail_numbers[key] = []
		nodes_count[key] = []
		flight_nums_inverted[key] = []
		crossings_inverted[key] = 0
		tail_numbers_inverted[key] = []
		node_inverted[key] = []
		for eq_type in equip_type:
			if str(eq_type) == str(key):
				ind = equip_type[equip_type == eq_type].index[0]
				engine_counts[key] = engine_count[ind]
				wake_turbs[key] = wake_turb[ind]
				type_name[key] = model[ind]
				engine_type_dict[key] = engine[ind]
				man_dict[key] = man[ind]
				break
			else:
				continue
		for line in infile_lines:
			data = line.split(',')  # Split the line by commas
			equip = data[10]  # Get the equipment type from the line
			if equip not in Equipment_jets:
				continue
			if equip == key:
				crossings[key] = crossings.get(key, 0) + 1
				flight_num = data[1]
				if flight_num not in flight_nums[key]:
					flight_nums[key].extend([flight_num])
				sta = data[9]
				date = data[0]
				flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_flights.csv', sep=",")
				flight = flight_data['flight_id']
				tailnumber = flight_data['aircraft_id']
				for i, f_id in enumerate(flight):
					if str(f_id) == str(flight_num):
						tailnum = tailnumber[i]
						break
					else:
						continue
				if tailnum not in tail_numbers[key]:
					tail_numbers[key].extend([tailnum])
				if sta not in nodes_count[key]:
					nodes_count[key].extend([sta])

		inv_file = 'output/inv_results/' + key + '_full_inv_results.txt'
		if not os.path.exists(inv_file):
			continue
		file = open('output/inv_results/' + key + '_full_inv_results.txt', 'r')
		for ll in file:
			l = ll.split(',')
			if str(l[-1]) == 'Forward Model':
				continue
			crossings_inverted[key] = crossings_inverted.get(key, 0) + 1
			f_num = l[1]
			if f_num not in flight_nums_inverted[key]:
				flight_nums_inverted[key].extend([f_num])
			sta = l[2]
			if sta not in node_inverted[key]:
				node_inverted[key].extend([sta])
			date = l[0]
			flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_flights.csv', sep=",")
			flight = flight_data['flight_id']
			tailnumber = flight_data['aircraft_id']
			for i, f_id in enumerate(flight):
				if str(f_id) == str(f_num):
					tailnum = tailnumber[i]
					break
				else:
					continue
			if tailnum not in tail_numbers_inverted[key]:
				tail_numbers_inverted[key].extend([tailnum])

	final_table = open('paper_table_jets.txt', 'w')
	# Write the header row
	final_table.write('Equipment,Manufacturer,Model,Engine Type,Engine Count,Wake Turbulence Category,Total Flights,Inverted Flights,Total Crossings,Inverted Crossings,Total Tail Numbers,Inverted Tail Numbers,Total Nodes,Inverted Nodes\n')
	for eq in Equipment_jets:

		final_table.write(eq +',' + str(man_dict[eq]) + ',' + str(type_name[eq]) + ',' + str(engine_type_dict[eq]) + ',' + str(engine_counts[eq]) + ',' + str(wake_turbs[eq]) + ',' + str(len(flight_nums[eq])) + ',' + str(len(flight_nums_inverted[eq])) + ',' + str(crossings[eq]) + ',' + str(crossings_inverted[eq]) + ',' + str(len(tail_numbers[eq])) + ',' + str(len(tail_numbers_inverted[eq])) + ',' + str(len(nodes_count[eq])) +  ',' + str(len(node_inverted[eq])) + '\n')

	final_table.close()

