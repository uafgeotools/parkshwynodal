import os

file = open('comb_output.txt', 'r')
list_list = []
lines = file.readlines()

for i,line in enumerate(lines):
	split_line = line.split(",")
	equip = split_line[-2]
	date = split_line[0]
	flight_num = split_line[1]
	time_stamp = split_line[3]
	sta = split_line[2]
	if equip[-1] == '_':
		equip = equip[:-1]
	list_list.append([equip, date, flight_num, time_stamp, sta])


input_dir  = '/scratch/irseppi/nodal_data/plane_info/inverse_final_database_NGT/'


for file_name in os.listdir(input_dir):
	if not file_name.endswith('.pdf'):
		continue
	split_name = file_name.split("_")
	equip = split_name[0]
	date = split_name[1]
	flight_num = split_name[2]
	time_stamp = split_name[3]
	sta = split_name[4]
	if [equip, date, flight_num, time_stamp, sta] in list_list:
		continue
	else:

		print(equip, date, flight_num, time_stamp, sta)

'''
main_dir = '/scratch/irseppi/nodal_data/plane_info/inverse_final_database_NGT/'
bad_dir = '/scratch/irseppi/nodal_data/plane_info/inverse_final_database_NGT/Bad_Fits/'

for fn in os.listdir(bad_dir):
	if not fn.endswith('.pdf'):
		continue
	main_path = os.path.join(main_dir, fn)
	print(main_path)
	if os.path.isfile(main_path):
		try:
			os.remove(main_path)
			print("Removed:", main_path)
		except OSError as e:
			print("Failed to remove", main_path, ":", e)
import os
def combine_text_files(input_files, output_file):

	with open(output_file, 'w') as outfile:
		for fname in input_files:
			equip = fname[0:4]
			with open(fname) as infile:
				for line in infile:
					if line.strip():
						outfile.write(line.rstrip('\n') + equip + ',\n')

input_dir = '/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results_ngt/'
output_file = 'comb_output.txt'
input_files = []

for file_name in os.listdir(input_dir):
    input_files.append(file_name)
combine_text_files(input_files, output_file)
'''
