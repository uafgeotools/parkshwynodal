#Read in the last two columns of a file  of the format 2019-02-11T22:46:34.996000Z|2019-03-09T22:46:34.994000Z
#Count the amount of days between the two dates and store the result in a list
#Average the list and print the result
from datetime import datetime, timezone
import numpy as np
import pandas as pd
file = pd.read_csv('input/parkshwy_nodes.txt', sep="|")

start_time = file['StartTime']
end_time = file['EndTime']

date_diffs = []

for i in range(len(start_time)):
    start_date =  datetime.strptime(start_time[i], "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)
    end_date = datetime.strptime(end_time[i], "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)

    diff = (end_date - start_date).days
    date_diffs.append(diff)
average_diff = np.mean(np.array(date_diffs))
print(f'Average difference in days: {average_diff}')
