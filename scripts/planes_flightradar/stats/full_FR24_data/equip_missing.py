import pandas as pd
import os
import numpy as np
from pathlib import Path

flight_files = []

for month in range(2,4):
	if month == 2:
		for day in range(11,29):
			day = str(day)
			# assign directory
			f = '/scratch/irseppi/nodal_data/flightradar24/201902'+month+day+'_flights.csv'

			# checking if it is a file
			if os.path.isfile(f):
				flight_files.append(f)
	elif month == 3:
		for day in range(1, 10):
			day = str(day)
			# assign directory
			f = '/scratch/irseppi/nodal_data/flightradar24/2019030'+day+'_flights.csv'
			
			# checking if it is a file
			if os.path.isfile(f):
				flight_files.append(f)

		for day in range(10, 27):
			day = str(day)
			# assign directory
			f = '/scratch/irseppi/nodal_data/flightradar24/201903'+day+'_flights.csv'
			# checking if it is a file
			if os.path.isfile(f):
				flight_files.append(f)

for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	e = flight_data['equip']
	fl = flight_data['aircraft_id']
	for line in range(len(e)):
		if e[line] == 'CNA':
			print(flight_file,fl[line],e[line]) 
		elif e[line] == 'BES':
			print(flight_file,fl[line],e[line])
		elif e[line] == 'DH8':
			print(flight_file,fl[line],e[line])
		elif e[line] == 'DH1':
			print(flight_file,fl[line],e[line])
		elif e[line] == 'CNC':
			print(flight_file,fl[line],e[line])
		elif e[line] == '73F':
			print(flight_file,fl[line],e[line])
		elif e[line] == 'D6F':
			print(flight_file,fl[line],e[line])
		elif e[line] == 'BE1':
			print(flight_file,fl[line],e[line])
		elif e[line] == 'KODI':
			print(flight_file,fl[line],e[line])


