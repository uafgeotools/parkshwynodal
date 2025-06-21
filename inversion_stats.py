import pandas as pd
import os

# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')
equip_count_dict = {}
tail_dict = {}
for li in file_in.readlines():
    text = li.split(',')
    date = text[0]
    flight_num = text[1]
    equip = text[10]
    sta = text[9]

    folder_spec = equip + '_spec_c'
    exist_dir = '/scratch/irseppi/nodal_data/plane_info/' + folder_spec +'/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'
    if os.path.exists(exist_dir):
    
        flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + date + '_flights.csv', sep=",")
        flight = flight_data['flight_id']
        flight = flight.values.tolist()
        tailnumber = flight_data['aircraft_id']

        for i,fly in enumerate(flight):
            if float(fly) == float(flight_num):
                tail = tailnumber[i]
                if tail not in tail_dict:
                    tail_dict[tail] = []
                    equip_count_dict[equip] = equip_count_dict.get(equip, 0) + 1
                break
    else:
        continue

print(equip_count_dict)