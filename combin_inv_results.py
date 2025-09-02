#read in every text file in a directory and combine all lines from all text files into a single text file
import os 
output = open('nodal_flight_inversion_results_ngt.txt', 'a')
base_dir = 'output/inv_results_no_g_truth/'
for file_name in os.listdir(base_dir):
	file_path = os.path.join(base_dir, file_name)
	if file_name.endswith('.txt'):
		with open(file_path, 'r') as f:
			lines = f.readlines()
			# Do something with the lines
			output.writelines(lines)
output.close()