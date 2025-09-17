import sys
import fileinput
import os
import shutil
import pandas as pd
from datetime import datetime

#############################################################################################################################

def modify_file(input_file_name, output_file_name):
	"""
	Takes all the letters in the input file and makes them uppercase, then writes this modified content to the output file.

	Args:
		input_file_name (str): The path of the input file.
		output_file_name (str): The path of the output file.

	Returns:
		None
	"""

	# Read the input file
	with open(input_file_name, 'r') as file:
		content = file.read()

	# Modify the content
	modified_content = content.upper()

	# Write the modified content to the output file
	with open(output_file_name, 'w') as file:
		file.write(modified_content)

####################################################################################################################

def station_subset(filename, steps, outputname):
	"""
	Subset stations from the input file based on the given steps and save the result to the output file.

	Parameters:
	- filename (str): The path of the input file.
	- steps (int): The amount of steps to subset stations by. (ie. every 4th station is in the subset)
	- outputname (str): The filename for the output file.
	"""

	output = open(outputname, 'w')
	with open(filename) as handle:
		for lineno, line in enumerate(handle):
			if lineno % steps == 0:
				output.write(line)
	output.close()

#######################################################################################


def replace(filename, old_string, new_string):
	"""
	Replace all occurrences of 'old_string' with 'new_string' in the given file.

	Args:
		filename (str): The path of the file to be modified.
		old_string (str): The string to be replaced.
		new_string (str): The string to replace the old_string with.
	"""

	for i, line in enumerate(fileinput.input(filename, inplace=1)):
		sys.stdout.write(line.replace(old_string, new_string))

	# Example usage:
	# replace('filename.site', '', ' "')
	# or replace('#', '\n#')

############################################################################################################################

def round_replace(filename, col_2round, precision, m2km):
	"""
	Replace the values, in meters, in a specific column of a text file with rounded values, in either meters or kilometers.

	Args:
		filename (str): The path of the text file.
		col_2round (int): The column index to round.
		precision (int): The number of digits to round to.
		m2km (bool): Determines whether to convert the rounded value to kilometers.
			if m2km = 0 - replace number with rounded number
			if m2km = 1 - replace number with rounded number converted to km
	
	Returns:
		None
	"""

	col_2round = int(col_2round)
	precision = int(precision)

	for i, line in enumerate(fileinput.input(filename, inplace=1)):
		val = line.split()

		if m2km == False:
			new_val = round(float(val[col_2round]), precision)
			sys.stdout.write(line.replace(str(val[col_2round]), str(new_val)))

		if m2km == True:
			new_val = round(float(val[col_2round]) / 1000, precision)
			sys.stdout.write(line.replace(str(val[col_2round]), str(new_val)))

#################################################################################################################################

def rename_file(flight_name):
	"""
	Renames all files in the specified flight collection so that the flight name is appended
	to the beginning of the file name, instead of a folder label.

	Args:
		flight_name (str): The name of the flight collection.

	Returns:
		None
	"""

	os.getcwd()
	collection = flight_name + '/'
	for i, filename in enumerate(os.listdir(collection)):
		for p, fil in enumerate(os.listdir(collection+filename)):
			os.rename(collection + filename + '/' + fil, collection + filename +'_'+ fil)

########################################################################################################
			
def extract_flight(equipment):
	"""
	Extracts all rows from the 'nodes_crossings_db_UTM.txt' file whith the designated equipment type 
	and prints them into an individual file labeled with the equipment type.
	
	Args:
		equipment (str): The equipment type to extract from the file.

	Returns:
		None
	"""

	input = open('input/nodes_crossings_db_UTM.txt','r')
	output = open('input/nodes_crossings_db_UTM_'+str(equipment)+'.txt','w')

	for line in input.readlines():
		val = line.split(',')
		if str(val[7][0:4]) == str(equipment):
			
			output.write(line)
		
	input.close()
	output.close()	

########################################################################################################
	
def extract_col(input_file, output_file, col, split_str):
	"""
	Extracts a specific column from a text file and writes it to another file.

	Args:
		input_file (str): The path to the input text file.
		output_file (str): The path to the output file where the extracted column will be written.
		col (int): The index of the column to extract (0-based index).
		split_str (str): The string that splits the text file into columns.

	Returns:
		None
	"""

	i = int(col)
	with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
		for line in f_in:
			line = line.split(split_str)
			print(line[i])
			f_out.write(line[i])

#########################################################################################

def date_round(input_file, output_file):
	"""
	Rounds the seconds of each timestamp in the input file to remove the milliseconds
	and writes the rounded timestamps to the output file.

	Args:
		input_file (str): The path to the input file containing timestamps.
		output_file (str): The path to the output file where rounded timestamps will be written.

	Returns:
		None
	"""

	# Remove the milliseconds from the timestamp
	with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
		for line in f_in:
			dt = datetime.strptime(line.strip(), '%Y-%m-%dT%H:%M:%S.%f')
			rounded_dt = dt.replace(second=round(dt.second))
			f_out.write(rounded_dt.strftime('%Y-%m-%d %H:%M:%S') + '\n')

########################################################################################

def count_flight(input_file, col_f, output_file, designator): 
	"""
	Counts the number of flights with a specific aircraft type designator in the input file.

	Args:
		input_file (str): Path to the input file.
		col_f (int): Column index where the flight equipment is located.
		output_file (str): Path to the output file.
		designator (str): Aircraft type designator to search for.

	Returns:
		None
	"""

	text = open(input_file, 'r')
	f = open(output_file, 'w')
	i = int(col_f)

	flight_data = pd.read_csv('20231010_Aircraft_UA_Fairbanks.csv', sep=",")
	eq = flight_data['TypeDesignator']
	des = flight_data['Description']

	count = 0
	for line in text.readlines():
		val = line.split(',')
		equip = val[i]
		for l in range(len(eq)):
			if str(eq[l]) == str(equip[0:4]) and str(des[l]) == designator:
				count = count + 1
				f.write(eq[l]+'\n')
	f.write(str(count))
	f.close()

#######################################################################################

def print_eq():
	"""
	Prints the aircraft type designator for each line in the 'nodes_crossings_db_UTM.txt' file.
	"""

	text = open('nodes_crossings_db_UTM.txt', 'r')

	for line in text.readlines():
		val = line.split(',')
		equip = val[6]
		print(equip)

#######################################################################################

def comb_lines(filename):
	"""
	Combines lines in the input file that are part of the same record.

	Args:
		filename (str): The path of the file to be modified.

	Returns:
		None
	"""

	with open(filename, "r") as file:
		lines = file.readlines()

	new_lines = []
	prev_line = ""

	for line in lines:
		if line.startswith(" "):
			prev_line += line.strip()
		else:
			new_lines.append(prev_line)
			prev_line = line.strip()

	# Append the last line
	new_lines.append(prev_line)

	# Write the modified lines back to the file
	with open(filename, "w") as file:
		file.write("\n".join(new_lines))

############################################################################################################

def order_rows_by_column(filename, col,split_symbol=',',option=1):
	"""
	Function to take one column of a text file and order it in increasing order. 
	Using that column, all rows will be rearranged in order.

	Args:
		filename (str): Path to the input file.
		col (int): The column index to sort by (0-based index).
		split_symbol (str): The symbol used to split the columns in the file. Default is ','.
		option (int): Determines the sorting method.
		1 - Sorts the values in the specified column and rearranges the lines based on the sorted values.
		2 - Sorts the entire file based on the third column.
	"""


	if option == 1:
		with open(filename, 'r') as file:
			lines = file.readlines()

		# Extract the values from the specified column
		values = []
		for line in lines:
			columns = line.split(split_symbol)

			values.append(columns[col].strip())
		print(values)
		# Sort the values based on the column
		sorted_values = sorted(values)

		# Rearrange the lines based on the sorted values
		rearranged_lines = []
		for value in sorted_values:
			for line in lines:
				if value in line:
					rearranged_lines.append(line)
					break

		# Write the rearranged lines back to the file
		with open(filename, 'w') as file:
			file.writelines(rearranged_lines)
			
	if option == 2: 
		# Read the file
		with open(filename, 'r') as file:
			lines = file.readlines()

		# Sort the lines based on the third column
		sorted_lines = sorted(lines, key=lambda line: line.split(',')[3])

		# Write the sorted data back to the file
		with open('sorted_'+filename, 'w') as file:
			file.writelines(sorted_lines)
############################################################################################################

def check_matching_values(file1, col1, file2, col2):
	"""
	Compares values in specific columns of two text files and prints the values if they differ.

	Args:
		file1 (str): Path to the first input file.
		col1 (int): Column index in the first file to compare (0-based index).
		file2 (str): Path to the second input file.
		col2 (int): Column index in the second file to compare (0-based index).
	"""

	with open(file1, 'r') as f1, open(file2, 'r') as f2:
		lines1 = f1.readlines()
		lines2 = f2.readlines()

	for i, (line1, line2) in enumerate(zip(lines1, lines2)):
		columns1 = line1.split(',')
		value1 = columns1[col1].strip()

		columns2 = line2.split(',')
		value2 = columns2[col2].strip()
		print(value1, value2)
		if value1 != value2:
			print(f"Row {i+1} in {file1} and {file2} have different values.")
############################################################################################################
		
def cojoin_columns(file1, start_col1, end_col1, file2, start_col2, end_col2, output_file):
	"""
	Function to take two text files, extract specific columns from each file, and cojoin them into a new file.

	Args:
		file1 (str): Path to the first input file.
		start_col1 (int): Starting column index for the first file.
		end_col1 (int): Ending column index for the first file.
		file2 (str): Path to the second input file.
		start_col2 (int): Starting column index for the second file.
		end_col2 (int): Ending column index for the second file.
		output_file (str): Path to the output file.

	Returns:
		None
	"""	

	with open(file1, 'r') as f1, open(file2, 'r') as f2, open(output_file, 'w') as output:
		lines1 = f1.readlines()
		lines2 = f2.readlines()


		for line1, line2 in zip(lines1, lines2):
			columns1 = line1.split(',')
			values1 = [columns1[i].strip() for i in range(start_col1, end_col1 + 1)]

			columns2 = line2.split(',')
			values2 = [columns2[i].strip() for i in range(start_col2, end_col2 + 1)]

			cojoined_line = ','.join(values1 + values2) + '\n'
			output.write(cojoined_line)

##############################################################################################################

def delete_empty_file(file_path):
    """
    Checks if a given file is empty and deletes it if it is.

    Args:
        file_path (str): The path to the file to be checked.
    """

    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' does not exist.")
        return

    try:
        if os.path.getsize(file_path) == 0:
            os.remove(file_path)
            print(f"File '{file_path}' was empty and has been deleted.")
        else:
            print(f"File '{file_path}' is not empty.")
    except OSError as e:
        print(f"Error processing file '{file_path}': {e}")

###############################################################################################################

def remove_files_with_no_picks(base_dir='/home/irseppi/REPOSITORIES/parkshwynodal/output/'):
	"""
	Removes files that do not contain any data picks from users.
	Args:
		base_dir (str): The base directory where the files are located. Default is '/home/irseppi/REPOSITORIES/parkshwynodal/output/'.
	"""

	# loop through the directories in the directory
	for dir_name in os.listdir(base_dir):
		dir_path = os.path.join(base_dir, dir_name)
		if os.path.isdir(dir_path):
			for equip_dir in os.listdir(dir_path):
				equip_dir_path = os.path.join(dir_path, equip_dir)
				if os.path.isdir(equip_dir_path):
					for date_dir in os.listdir(equip_dir_path):
						date_dir_path = os.path.join(equip_dir_path, date_dir)
						if os.path.isdir(date_dir_path):
							for flight_dir in os.listdir(date_dir_path):
								flight_dir_path = os.path.join(date_dir_path, flight_dir)
								if os.path.isdir(flight_dir_path):
									for sta_dir in os.listdir(flight_dir_path):
										sta_dir_path = os.path.join(flight_dir_path, sta_dir)
										if os.path.isdir(sta_dir_path):
											for file_name in os.listdir(sta_dir_path):
												file_path = os.path.join(sta_dir_path, file_name)
												delete_empty_file(file_path)

################################################################################################################

def remove_dir_with_no_picks(base_dir='/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/'):
	"""
	Removes directories that do not contain any data picks from users.
	Args:
		base_dir (str): The base directory where the files are located. Default is '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/'.
	"""
	for dir_name in os.listdir(base_dir):
		dir_path = os.path.join(base_dir, dir_name)
		# Walk through the directory tree from the bottom up
		for root, dirs, files in os.walk(dir_path, topdown=False):
			# If the directory is empty (no files and no subdirectories), delete it
			if not files and not dirs:
				os.rmdir(root)
				print(f"Deleted empty directory: {root}")

################################################################################################################

def clean_inv_results():
	"""
	Cleans up the formatting of files in the 'output/inv_results/' directory by replacing specific patterns.
	Replaces all instances of '\n ' with ' ', '  ' with ' ', '[ ' with '[', and ' ]' with ']'.
	"""

	# specify the directory containing the files to be cleaned
	dir_path = '/home/irseppi/REPOSITORIES/parkshwynodal/output/inv_results_no_g_truth_320/'

	# loop through the directories in the directory
	for dir_name in os.listdir(dir_path):
		filename = os.path.join(dir_path, dir_name)

		for line in fileinput.input(filename, inplace=1):
			line = line.replace('\n ', ' ') #.replace('  ', ' ').replace('[ ', '[').replace(' ]', ']')
			sys.stdout.write(line)

################################################################################################################

def combine_text_files(input_files, output_file):
	"""
	Combines multiple text files into a single output file.

	Args:
		input_files (list): List of input file names to combine.
		output_file (str): The path to the output file.

	Example usage:
		files_to_combine = []
		for file in os.listdir('output/inv_results/'):
			if file.endswith('.txt'):
				files_to_combine.append(file)
		output_filename = 'combined_python.txt'
		combine_text_files(files_to_combine, output_filename)
	"""
	with open(output_file, 'w') as outfile:
		for fname in input_files:
			equip = fname[0:4]
			with open(os.path.join('output/inv_results/', fname)) as infile:
				for line in infile:
					if line.strip():
						outfile.write(line.rstrip('\n') + equip + ',\n')

#############################################################################################################################

def change_directory_structure(BASE_DIR, BASE_DIR2):
	"""
	Changes the directory structure by copying files from the source directory to the target directory while maintaining the
	original structure.
	Args:
		BASE_DIR (str): The source base directory.
		BASE_DIR2 (str): The target base directory.
	"""
	for file_name in os.listdir(BASE_DIR):
		sig_file = os.path.join(BASE_DIR, file_name)
		sig_hold = os.path.join(BASE_DIR2, file_name)
		for eq in os.listdir(sig_file):
			eq_file = os.path.join(sig_file, eq)
			for date_file in os.listdir(eq_file):
				full_path = os.path.join(eq_file, date_file)
				for flight_file in os.listdir(full_path):
					flight_path = os.path.join(full_path, flight_file)
					for sta_file in os.listdir(flight_path):
						sta_path = os.path.join(flight_path, sta_file)
						for file in os.listdir(sta_path):
							image_file = os.path.join(sta_path, file)
							new_image_file = os.path.join(sig_hold,file)
							os.makedirs(os.path.dirname(new_image_file), exist_ok=True)
							shutil.copy2(image_file, new_image_file)

#############################################################################################################################

def combine_all_text_files_in_dir(input_dir, output_file):
	"""
	Reads all text files in a directory and combines all lines from all text files into a single output text file.

	Args:
		input_dir (str): The directory containing the text files.
		output_file (str): The path to the output file.
	"""
	with open(output_file, 'a') as output:
		for file_name in os.listdir(input_dir):
			file_path = os.path.join(input_dir, file_name)
			if file_name.endswith('.txt'):
				with open(file_path, 'r') as f:
					lines = f.readlines()
					output.writelines(lines)

#############################################################################################################################

def compare_directories(dir1, dir2, dir3):
	"""
	Compares files in three directories and prints out any files that are missing in either of the
	second or third directories compared to the first directory.
	Args:
		dir1 (str): The path of the first directory.
		dir2 (str): The path of the second directory.
		dir3 (str): The path of the third directory.
	Returns:
		int: The total number of files checked across all three directories.

	EXAMPLE USAGE:
	AIRCRAFT_TYPES = 0
	for dirc in os.listdir('input/Data_Picks'):
		dirc_list = []
		new_dirc = os.path.join('input/Data_Picks', dirc)
		for dirc in os.listdir(new_dirc):
			dirc_list.append(os.path.join(new_dirc, dirc))
		AIRCRAFT_TYPES += compare_directories(dirc_list[0], dirc_list[1], dirc_list[2])
	print(AIRCRAFT_TYPES)
	"""
	
	FILE_COUNT = 0
	for root, _, files in os.walk(dir1):
		for filename in files:
			print(root)
			file1 = os.path.join(root, filename)
			file2 = os.path.join(dir2, os.path.relpath(file1, dir1))
			file3 = os.path.join(dir3, os.path.relpath(file1, dir1))
			FILE_COUNT += 3
			if not os.path.exists(file2):
				print(f"Missing in {dir2}: {file2}")
			if not os.path.exists(file3):
				print(f"Missing in {dir3}: {file3}")
	return FILE_COUNT

#############################################################################################################################
