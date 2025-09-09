import pandas as pd
import matplotlib.pyplot as plt
import datetime
import imageio
from obspy.geodetics import gps2dist_azimuth

def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

min_lon = -150.5
max_lon = -148.5
min_lat = 62.227
max_lat = 64.6


# Load the seismometer location data
seismo_data = pd.read_csv('nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_stations = seismo_data['Latitude']
sta = seismo_data['Station']

	
f = '/scratch/irseppi/nodal_data/flightradar24/20190314_positions/20190314_533595343.csv'


flight_data = pd.read_csv(f, sep=",")
flight_latitudes = flight_data['latitude']
flight_longitudes = flight_data['longitude']
time = flight_data['snapshot_id']
speed = flight_data['speed']
alt = flight_data['altitude']

flight_num = '533595343'	

flight_data['timestamp'] = pd.to_datetime(time)

# Sort flight data by timestamp
flight_data = flight_data.sort_values('timestamp')

fig = plt.figure(figsize=(24,30))
# Create a scatter plot of the seismometer locations
plt.scatter(seismo_longitudes, seismo_latitudes, color='red')

# List to store each frame of the GIF
images = []

for i in range(len(flight_data)):
	for sd in range(len(seismo_data)):	
		dist = distance(seismo_latitudes[sd], seismo_longitudes[sd], flight_latitudes[i], flight_longitudes[i])
		if dist <= 2:
			# Add the flight timestamps to the plot
			plt.scatter(flight_data.iloc[i]['longitude'], flight_data.iloc[i]['latitude'], color='blue')
			#Label time stamps with time
			plt.text(flight_longitudes[i], flight_latitudes[i], flight_data['timestamp'][i])
			# Set labels and title
			plt.xlim(min_lon, max_lon)
			plt.ylim(min_lat, max_lat)
			plt.xlabel('Longitude')
			plt.ylabel('Latitude')
			plt.title('Flight: 533595343')
			# Save the plot as a PNG file
			plt.savefig(f'frame_{i}.png')

			# Add the PNG file to the list of images
			images.append(imageio.imread(f'frame_{i}.png'))

			
		else:
			continue
# Create a GIF from the images with 5 frames per second (one frame every 0.2 seconds)
imageio.mimsave('flight.gif', images, fps=4)
