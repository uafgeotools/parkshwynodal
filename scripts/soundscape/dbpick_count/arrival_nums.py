from datetime import datetime
from collections import Counter
import matplotlib.pyplot as plt

# Read the timestamps from the text file
with open('/REPOSITORIES/parkshwynodal/input/dates.txt', 'r') as f:
    timestamps = f.read().splitlines()

# Convert the timestamps to datetime objects
datetimes = [datetime.strptime(ts, '%Y-%m-%d %H:%M:%S') for ts in timestamps]

# Filter the datetimes to only include those between 8 and 15 UTC
filtered_datetimes = [dt for dt in datetimes if 8 >= dt.hour or dt.hour > 15]

# Count the number of dates for each day
date_counts = Counter(dt.date() for dt in filtered_datetimes)

# Create a histogram from the date counts
plt.bar(date_counts.keys(), date_counts.values())
plt.title('Day Time Arrival Count')
plt.xlabel('Date')
plt.ylabel('Count')
plt.show()
	

# Read the timestamps from the text file
with open('/REPOSITORIES/parkshwynodal/input/dates.txt', 'r') as f:
    timestamps = f.read().splitlines()

# Convert the timestamps to datetime objects
datetimes = [datetime.strptime(ts, '%Y-%m-%d %H:%M:%S') for ts in timestamps]

# Filter the datetimes to only include those between 8 and 15 UTC
filtered_datetimes = [dt for dt in datetimes if 8 <= dt.hour < 15]

# Count the number of dates for each day
date_counts = Counter(dt.date() for dt in filtered_datetimes)

# Create a histogram from the date counts
plt.bar(date_counts.keys(), date_counts.values())
plt.title('Night Time Arrival Count')
plt.xlabel('Date')
plt.ylabel('Count')
plt.show()

plt.hist(datetimes, bins=43)
plt.title('All Time Arrival Count')
plt.xlabel('Date')
plt.ylabel('Count')
plt.show()
