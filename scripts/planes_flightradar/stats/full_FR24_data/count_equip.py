import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np

# Get a list of all the csv files in the 'YYYYDDMM_positions' folder
files = glob.glob('/scratch/irseppi/nodal_data/flightradar24/*_flights.csv')

# Initialize an empty DataFrame
dfs = []

# Loop through the files and append the data to the DataFrame
for file in files:
    data = pd.read_csv(file)
    dfs.append(data)

# Concatenate all the dataframes in the list
df = pd.concat(dfs)
print(len(df))
# Count the occurrences of each type per day
df_counts = df['equip']
# Count the number of different equipment types
num_equipment_types = df_counts.nunique()
print(num_equipment_types)

#How many flights dont have equipment type
print(df_counts.isnull().sum())

file_2 = 'all_station_crossing_db_updated111.txt'
# Read in file 2 and count all of the flight numbers column 2 but a single occurrence of each flight number
file2 = pd.read_csv(file_2, sep=",")
flight = file2['527831677']
print(flight.nunique())

date = file2['20190211']  # Assuming the first column is the day column
unique_flights_per_day = flight.groupby(date).nunique()
print(unique_flights_per_day)
print(np.mean(unique_flights_per_day))
print(np.mean(unique_flights_per_day)-np.min(unique_flights_per_day))
print(np.max(unique_flights_per_day)-np.mean(unique_flights_per_day))

