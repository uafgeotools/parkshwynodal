import os
from obspy import read

file_in = open('input/parkshwy_nodes.txt', 'r+')
folder = '/scratch/irseppi/500sps/'
file_out = open('input/parkshwy_nodes_corrected.txt', 'w')
# Read all lines from the file once
file_lines = file_in.readlines()

# Preprocess the folder structure to avoid repeated os.listdir calls
start_end = {}
for date_f in os.listdir(folder):
    date_path = os.path.join(folder, date_f)
    for station_f in os.listdir(date_path):
        if station_f.endswith('DPZ.msd'):
            sta = station_f.split('_')[1]                    
            month = date_f.split('_')[1]
            day = date_f.split('_')[2]
            da = [month, day]
            num = int(''.join(da))
            if sta not in start_end:
                start_end[sta] = []
            start_end[sta].extend([num])

file_structure = {}
for sta, date in start_end.items():
    min_day = str(min(date))
    max_day = str(max(date))
    date_1 = '2019_0' + str(min_day[0]) + '_' + str(min_day[1:])
    date_2 = '2019_0' + str(max_day[0]) + '_' + str(max_day[1:])
    file_structure[sta] = [os.path.join(folder, date_1, 'ZE_' + sta + '_DPZ.msd'), os.path.join(folder, date_2, 'ZE_' + sta + '_DPZ.msd')]
x = 0

for sta in file_structure.keys():
    path = file_structure[sta]
    st = read(path[0])
    start_time = st[0].stats.starttime
    st = read(path[1])
    end_time = st[0].stats.endtime
    for line in file_lines:
        text = line.split('|')
        if text[1] == sta:
            x += 1
            print(sta, x)
            text[6] = str(start_time)
            text[7] = str(end_time)
            updated_line = '|'.join(text)

            # Update the line in file_lines
            file_lines[file_lines.index(line)] = updated_line + '\n'
            break
        else:
            continue
file_out.writelines(file_lines)
file_in.close()
file_out.close()
