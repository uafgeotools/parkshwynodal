# Import the necessary modules
import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

# Define the parameters
#Using PH5ws for 02/13, 02/20, 02/23, and 03/01
network = "ZE"
day = "2019-02-13"
#day = "2019-02-20"
#day = "2019-02-23"
#day = "2019-03-01"
report_number = "19-007"
client = Client("http://service.iris.edu/ph5ws/", service_mappings={
    "dataselect": "http://service.iris.edu/ph5ws/dataselect/1",
    "station": "http://service.iris.edu/ph5ws/station/1"
})

# Get the inventory of the network
inventory = client.get_stations(network=network, level="response")

# Loop through the stations and channels
for station in inventory[0]:
	for channel in station:
		# Define the start and end times of the data
		start_time = UTCDateTime(day + "T00:00:00")
		end_time = UTCDateTime(day + "T23:59:59")

		segment_start = start_time
		segment_end = end_time

		# Try to download the data for each segment
		try:
			waveform = client.get_waveforms(network, station.code, "*", channel.code, segment_start, segment_end)

			# Save the waveform as a miniseed file with the station and component as a suffix
			filename = "ZE_"+station.code + "_" + channel.code + ".msd"
			waveform.write(filename, format="MSEED")

			# Print a message to indicate success
			print(f"Downloaded {filename}")
		except:
			print('not there')
