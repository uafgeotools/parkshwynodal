import pandas as pd

repo_path = '/home/irseppi/REPOSITORIES/parkshwynodal/'
flightradar24_data_path = '/scratch/irseppi/nodal_data/flightradar24/'

jets = False
with open(repo_path + 'input/node_crossings_db_UTM.txt', 'r') as infile:
	infile_lines = infile.readlines()
if not jets:	
	Equipment = ['DH8A', 'B190','BE20','PC12','DH3T','C208','AT73','SW4','C441','B18T','B350','BE10','AS50','R44','C185','PA31','DHC2','GA8','C180','C182','C206','C172','PA32','PA46','CH7B','PA30','C46','BE35','PA18','PA34']
else:
	Equipment = ['B737','B738','B739','B77W','B772','B789','B788','B733','B763','A359','B77L','B744','E75S','B732','A332','B748','nan','CRJ2']

equip_data = pd.read_csv(repo_path + 'input/20231010_Aircraft_UA_Fairbanks.csv', sep=",")
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
		if str(eq_type) == str(key) and not str(eq_type) == 'nan':
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
		if equip == key :
			crossings[key] = crossings.get(key, 0) + 1
			flight_num = data[1]
			if flight_num not in flight_nums[key]:
				flight_nums[key].extend([flight_num])
			sta = data[9]
			date = data[0]
			flight_data = pd.read_csv(flightradar24_data_path + str(date) + '_flights.csv', sep=",")
			flight = flight_data['flight_id']
			tailnumber = flight_data['aircraft_id']
			for i, f_id in enumerate(flight):
				if str(f_id) == str(flight_num) and not str(equip) == 'Unknown':
					tailnum = tailnumber[i]
					break
				elif str(equip) == 'Unknown':
					tailnum = 'Unknown'
					break
				else:
					continue
			if tailnum not in tail_numbers[key]:
				tail_numbers[key].extend([tailnum])
			if sta not in nodes_count[key]:
				nodes_count[key].extend([sta])

	file = open(repo_path + 'output/NGT_flight_param_inv_DB.txt', 'r')
	for i, ll in enumerate(file):
		if i == 0:
			continue
		l = ll.split(',')

		if str(l[-2]) == 'Forward Model':
			continue
		if str(l[0]) != str(key):
			continue	

		crossings_inverted[key] = crossings_inverted.get(key, 0) + 1

		f_num = l[2]
		if f_num not in flight_nums_inverted[key]:
			flight_nums_inverted[key].extend([f_num])
		sta = l[3]
		if sta not in node_inverted[key]:
			node_inverted[key].extend([sta])
		date = l[1]
		flight_data = pd.read_csv(flightradar24_data_path + str(date) + '_flights.csv', sep=",")
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
	file.close()
if not jets:
	final_table = open(repo_path + 'paper_table.txt', 'w')
else:
	final_table = open(repo_path + 'paper_table_jets.txt', 'w')

# Write the header row
final_table.write('Equipment,Manufacturer,Model,Engine Type,Engine Count,Wake Turbulence Category,Total Flights,Inverted Flights,Total Crossings,Inverted Crossings,Total Tail Numbers,Inverted Tail Numbers,Total Nodes,Inverted Nodes\n')
for eq in Equipment:
	if str(eq) == 'nan':

		final_table.write(str('Unkown') +',' + str('Unkown') + ',' + str('Unkown') + ',' + str('Unkown') + ',' + str('Unkown') + ',' + str('Unkown') + ',' + str(len(flight_nums[eq])) + ',' + str(len(flight_nums_inverted[eq])) + ',' + str(crossings[eq]) + ',' + str(crossings_inverted[eq]) + ',' + str('Unkown') + ',' + str('Unkown') + ',' + str(len(nodes_count[eq])) +  ',' + str(len(node_inverted[eq])) + '\n')
	else:
		final_table.write(eq +',' + str(man_dict[eq]) + ',' + str(type_name[eq]) + ',' + str(engine_type_dict[eq]) + ',' + str(engine_counts[eq]) + ',' + str(wake_turbs[eq]) + ',' + str(len(flight_nums[eq])) + ',' + str(len(flight_nums_inverted[eq])) + ',' + str(crossings[eq]) + ',' + str(crossings_inverted[eq]) + ',' + str(len(tail_numbers[eq])) + ',' + str(len(tail_numbers_inverted[eq])) + ',' + str(len(nodes_count[eq])) +  ',' + str(len(node_inverted[eq])) + '\n')
#Make the last row totals of all the columns
total_flights = 0
total_inverted_flights = 0
total_crossings = 0
total_inverted_crossings = 0
total_tail_numbers = 0
total_inverted_tail_numbers = 0
total_nodes = 0
total_inverted_nodes = 0
for eq in Equipment:	
	total_flights += len(flight_nums[eq])
	total_inverted_flights += len(flight_nums_inverted[eq])
	total_crossings += crossings[eq]
	total_inverted_crossings += crossings_inverted[eq]
	total_tail_numbers += len(tail_numbers[eq])
	total_inverted_tail_numbers += len(tail_numbers_inverted[eq])
	total_nodes += len(nodes_count[eq])
	total_inverted_nodes += len(node_inverted[eq])
final_table.write('Totals,' + ',' + ',' + ',' + ',' + ',' + str(total_flights) + ',' + str(total_inverted_flights) + ',' + str(total_crossings) + ',' + str(total_inverted_crossings) + ',' + str(total_tail_numbers) + ',' + str(total_inverted_tail_numbers) + ',' + str(total_nodes) +  ',' + str(total_inverted_nodes) + '\n')
final_table.close()

#In dataset not only inverted
infile = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', 'r')
Equipment = ['DH8A', 'B190','BE20','PC12','DH3T','C208','AT73','SW4','C441','B18T','B350','BE10','AS50','R44','C185','PA31','DHC2','GA8','C180','C182','C206','C172','PA32','PA46','CH7B','PA30','C46','BE35','PA18','PA34']
Equipment_jet = ['B737','B738','B739','B77W','B772','B789','B788','B733','B763','A359','B77L','B744','E75S','B732','A332','B748','CRJ2']
total_nan_crossings = 0
prop_crossings = 0
total_jet_crossings = 0
index = 0
index_active = 0
for line in infile.readlines():
	data = line.split(',')
	equip = data[10]
	if equip not in Equipment and equip not in Equipment_jet and str(equip) != 'nan':
		print(equip)
	index += 1
	if str(equip) == 'nan':
		equip = 'Unknown'
		total_nan_crossings += 1
	if equip  in Equipment:
		prop_crossings += 1
	elif equip in Equipment_jet:
		total_jet_crossings += 1
print('Crossings in dataset (not only inverted):')
print('Total NaN Crossings: ' + str(total_nan_crossings))
print('Total Prop Crossings: ' + str(prop_crossings))
print('Total Jet Crossings: ' + str(total_jet_crossings))
print('Total Crossings: ' + str(total_nan_crossings + prop_crossings + total_jet_crossings))

#Inverted
infile = open(repo_path + 'output/NGT_flight_param_inv_DB.txt', 'r')
Equipment = ['DH8A', 'B190','BE20','PC12','DH3T','C208','AT73','SW4','C441','B18T','B350','BE10','AS50','R44','C185','PA31','DHC2','GA8','C180','C182','C206','C172','PA32','PA46','CH7B','PA30','C46','BE35','PA18','PA34']
Equipment_jet = ['B737','B738','B739','B77W','B772','B789','B788','B733','B763','A359','B77L','B744','E75S','B732','A332','B748','CRJ2']
total_nan_crossings = 0
prop_crossings = 0
total_jet_crossings = 0
index = 0
index_active = 0
for line in infile.readlines():
	data = line.split(',')
	equip = data[0]
	if equip not in Equipment and equip not in Equipment_jet and str(equip) != 'nan':
		print(equip)
	index += 1
	if str(equip) == 'nan':
		equip = 'Unknown'
		total_nan_crossings += 1
	if equip  in Equipment:
		prop_crossings += 1
	elif equip in Equipment_jet:
		total_jet_crossings += 1
print('Crossings in dataset that were inverted for:')
print('Total NaN Crossings: ' + str(total_nan_crossings))
print('Total Prop Crossings: ' + str(prop_crossings))
print('Total Jet Crossings: ' + str(total_jet_crossings))
print('Total Crossings: ' + str(total_nan_crossings + prop_crossings + total_jet_crossings))