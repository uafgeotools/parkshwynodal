import pandas as pd
import matplotlib.pyplot as plt
import glob
import matplotlib
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
import pytz
import obspy

flight_names=[]
dates = []

month = '02'
for day in range(11,29):
	day = str(day)
	# assign directory
	directory = '/scratch/irseppi/nodal_data/Plane_info/Plane_map/2019-' + month + '-' + day

	# iterate over files in directory
	for filename in os.listdir(directory):
		flight_name = filename[13:22]
		date = filename[4:12]

		flight_names.append(flight_name)
		dates.append(date)

month = '03'
for day in range(1, 10):
	day = '0' + str(day)
	# assign directory
	directory = '/scratch/irseppi/nodal_data/Plane_info/Plane_map/2019-' + month + '-' + day

	# iterate over files in directory
	for filename in os.listdir(directory):
		flight_name = filename[13:22]
		date = filename[4:12]

		flight_names.append(flight_name)
		dates.append(date)

for day in range(10, 27):
	day = str(day)
	# assign directory
	directory = '/scratch/irseppi/nodal_data/Plane_info/Plane_map/2019-' + month + '-' + day
	
	# iterate over files in directory
	for filename in os.listdir(directory):
		flight_name = filename[13:22]
		date = filename[4:12]

		flight_names.append(flight_name)
		dates.append(date)
				
f = open('designator_missing.txt','w')

for i in range(len(flight_names)):
	text = open('/scratch/irseppi/nodal_data/flightradar24/'+ dates[i] +'_flights.csv', "r")
	for line in text.readlines():
		val = line.split(',')
		if val[0] == flight_names[i] and val[1] == '0':
			val.append(dates[i])
			f.write(val[10]+','+val[0]+','+val[1]+','+val[2]+','+val[3]+','+val[4]+','+val[5]+','+val[6]+','+val[7]+','+val[8]+'\n')

f.close()				
#for equip with no designator
text = open('output/designator_missing.txt','r')
f = open('output/equp_missing.txt','w')
for line in text.readlines():
	val = line.split(',')
	print(val[4])
	if val[4] == '':
		f.write(val[0]+','+val[1]+','+val[2]+','+val[3]+','+val[4]+','+val[5]+','+val[6]+','+val[7]+','+val[8]+','+val[9])
f.close()
