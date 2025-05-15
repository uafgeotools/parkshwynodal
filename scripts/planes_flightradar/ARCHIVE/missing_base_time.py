import pandas as pd
import os
import numpy as np

file_p = 'output/C185_data_picks/inversepicks'
for date in os.listdir(file_p):
	f = os.path.join(file_p, date)
	for flight in os.listdir(f):
		flight_f = os.path.join(f, flight)
		for station in os.listdir(flight_f):
			sta = os.path.join(flight_f, station)
			for file in os.listdir(sta):
				file_path = os.path.join(sta, file)
				try:
					data = pd.read_csv(file_path, sep=",", header=None)
					column_3 = data.iloc[:, 2]
					for number in column_3:

						if str(number) == str(np.nan):
							print(file_path)
							break

				except:
					continue
