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

def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

def dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes):
	f = False
	for s in range(len(flight_latitudes)):
		for l in range(len(seismo_latitudes)):
			dist = distance(seismo_latitudes[l], seismo_longitudes[l], flight_latitudes[s], flight_longitudes[s])
			if dist <= 5:
				f = True
				break
			else:
				continue
	return f


min_lon = -150.5
max_lon = -148.5
min_lat = 62.227
max_lat = 64.6


flight_files=[]
filenames = []

# Load the seismometer location data
seismo_data = pd.read_csv('nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_stations = seismo_data['Latitude']
sta = seismo_data['Station']

for month in (2,4):
	if month == 2:
		month = '02'
		for day in range(11,29):
			day = str(day)
			# assign directory
			directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'

			# iterate over files in directory
			for filename in os.listdir(directory):
				filenames.append(filename)
				f = os.path.join(directory, filename)
				
				# checking if it is a file
				if os.path.isfile(f):
					flight_files.append(f)
	elif month == 3:
		month = '03'
		for day in range(1, 27):
			if day < 10:
				day = '0' + str(day)
				# assign directory
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'
			
				# iterate over files in directory
				for filename in os.listdir(directory):
					filenames.append(filename)
					f = os.path.join(directory, filename)
					
					# checking if it is a file
					if os.path.isfile(f):
						flight_files.append(f)
			else:
				day = str(day)
				# assign directory
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'
				
				# iterate over files in directory
				for filename in os.listdir(directory):
					filenames.append(filename)
					f = os.path.join(directory, filename)
					
					# checking if it is a file
					if os.path.isfile(f):
						flight_files.append(f)
				
# Initialize an empty list to store DataFrames
dfs = []

for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	fname = filenames[i]	
	flight_num = fname[9:18]
	date = fname[0:8]
	con = dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes)
	if con == True:	

		files = glob.glob(date+'_flights.csv')

		# Loop through the files and append the data to the DataFrame
		for file in files:
		    data = pd.read_csv(file)
		    dfs.append(data)

# Concatenate all the dataframes in the list
df = pd.concat(dfs)

# Extract the 'equipment type' column and count the occurrences of each type
equipment_counts = df['equip'].value_counts()

# Calculate the percentage of each equipment type
equipment_percent = equipment_counts / equipment_counts.sum() * 100

# Group equipment types that occur less than 0.5% of the time into 'Other'
#other_count = equipment_counts[equipment_percent < 0.5].sum()
#equipment_counts = equipment_counts[equipment_percent >= 0.5]

# Ensure 'Other' category always appears in the pie chart
#if 'Other' not in equipment_counts.index:
#equipment_counts = pd.concat([equipment_counts, pd.Series([0], index=['Other'])])

# Add the count of 'Other' equipment types
#equipment_counts['Other'] += other_count

# Plot the data as a bar graph
equipment_counts.plot(kind='bar')
plt.xlabel('Equipment Type')
plt.ylabel('Count')
plt.title('Occurrences of Equipment Types')
plt.show()

# Plot the data as a pie chart
#equipment_counts.plot(kind='pie', autopct='%1.1f%%')
#plt.title('Occurrences of Equipment Types')
#plt.show()
