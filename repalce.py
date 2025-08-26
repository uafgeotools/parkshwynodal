# go through each file in output/inv_results and replace all instances of '\n ' with ' ', repalce all '  ' with ' ', and repalce all '[ ' with '[' and all ' ]' with ']'
import sys
import os
import fileinput

dir_path = 'output/inv_results/'

# loop through the directories in the directory
for dir_name in os.listdir(dir_path):
    filename = os.path.join(dir_path, dir_name)

    for line in fileinput.input(filename, inplace=1):
        line = line.replace('\n ', ' ').replace('  ', ' ').replace('[ ', '[').replace(' ]', ']')
        sys.stdout.write(line)

	
