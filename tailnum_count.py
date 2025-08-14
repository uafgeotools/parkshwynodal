import pandas as pd

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

tail_num_dict = []
count = 0
for text in file_in.readlines():
    lines = text.split(',')
    date = lines[0]
    flight_num = lines[1]

    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + date + '_flights.csv', sep=",")
    flight = flight_data['flight_id']
    tailnumber = flight_data['aircraft_id']

    for i,fly in enumerate(flight):
        if float(fly) == float(flight_num):
            if tailnumber[i] not in tail_num_dict:
                tail_num_dict.append(tailnumber[i])
                count += 1

print(count)