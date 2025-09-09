import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime
import datetime
from datetime import datetime, timedelta
import pytz

# open file in read mode
text = open("nodal_ZE.arrival", "r")

day_data=[]	
night_data=[]
num = 0

for line in text.readlines():
	val = line.split()
	timestamps = UTCDateTime(float(val[1]))

	for x in range(11,29):

		if x <= 28:
			# Set the start and end times for the day and night periods
			day_start = UTCDateTime("2019-02-"+str(x-1)+"T15:00:00")
			day_end = UTCDateTime("2019-02-"+str(x)+"T08:08:00")
			night_start = UTCDateTime("2019-02-"+str(x)+"T08:00:00")
			night_end = UTCDateTime("2019-02-"+str(x)+"T15:00:00")

			# Filter the data to only include the day and night periods
			if day_start <= timestamps < day_end:
				day_data.append(timestamps)
			if night_start <= timestamps < night_end:
				night_data.append(timestamps)

	for y in range(2,28):
		day_start = UTCDateTime("2019-03-"+str(y-1)+"T15:00:00")
		day_end = UTCDateTime("2019-03-"+str(y)+"T08:00:00")
		night_start = UTCDateTime("2023-03-"+str(y)+"T08:00:00")
		night_end = UTCDateTime("2023-03-"+str(y)+"T15:00:00")

		# Filter the data to only include the day and night periods
		if day_start <= timestamps < day_end:
			day_data.append(timestamps)
		if night_start <= timestamps < night_end:
			night_data.append(timestamps)

	day_start = UTCDateTime("2019-02-28T15:00:00")
	day_end = UTCDateTime("2019-02-01T08:00:00")
	night_start = UTCDateTime("2023-03-01T08:00:00")
	night_end = UTCDateTime("2023-03-01T15:15:00")
	if day_start <= timestamps < day_end:
		day_data.append(timestamps)
	if night_start <= timestamps < night_end:
		night_data.append(timestamps)	

	num = num + 1
	print((num/2132047)*100)

# Calculate the arrival times for the day and night periods
day_arrivals = len(day_data)
night_arrivals = len(night_data)

# Compare the arrival times between day and night
print("# of Day arrivals:", day_arrivals)
print("# of Night arrivals:", night_arrivals)
print("Total # of picks", num)

# Convert UTCDateTime objects to datetime objects
day_dates = [datetime.utcfromtimestamp(date.timestamp) for date in day_data]
night_dates = [datetime.utcfromtimestamp(dates.timestamp) for dates in night_data]

# Define the local time zone
local_tz = pytz.timezone('US/Alaska')

# Convert UTC datetimes to local datetimes
local_day = [local_tz.normalize(day_date.replace(tzinfo=pytz.utc).astimezone(local_tz)) for day_date in day_dates]
local_night = [local_tz.normalize(night_date.replace(tzinfo=pytz.utc).astimezone(local_tz)) for night_date in night_dates]

# Create histogram
plt.hist(local_day, bins=43)
plt.title('Day Time Arrival Count')
plt.xlabel('Date')
plt.ylabel('Count')
plt.show()


plt.hist(local_night, bins=43)
plt.title('Night Time Arrival Count')
plt.xlabel('Date')
plt.ylabel('Count')
plt.show()
