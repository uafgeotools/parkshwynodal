import sys
import fileinput
import os
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

def order_rows_by_column(filename, col,split_symbol=','):
	"""
	Function to take one column of a text file and order it in increasing order. 
	Using that column, all rows will be rearranged in order.

	Args:
		filename (str): Path to the input file.
		col (int): The column index to sort by (0-based index).
		split_symbol (str): The symbol used to split the columns in the file. Default is ','.
	"""

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
