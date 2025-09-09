import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Get a list of all the csv files in the 'YYYYDDMM_positions' folder
files = glob.glob('*_flights.csv')
#files = [item for sublist in files for item in sublist]  # Flatten the list

# Initialize an empty DataFrame
dfs = []

# Loop through the files and append the data to the DataFrame
for file in files:
    data = pd.read_csv(file)
    # Add a 'date' column to the data
    date = os.path.basename(file).split('_')[0]
    data['date'] = pd.to_datetime(date, format='%Y%m%d').strftime('%m/%d')
    dfs.append(data)

# Concatenate all the dataframes in the list
df = pd.concat(dfs)

# Drop duplicate flights per day
#df = df.drop_duplicates(subset=['date', 'flight_number'])

# Count the occurrences of each type per day
df_counts = df.groupby(['date', 'equip']).size().unstack(fill_value=0)

# Calculate the percentage of each equipment type
equipment_percent = df_counts / df_counts.sum() * 100

# Group equipment types that occur less than 2% of the time into 'Other'
other_count = df_counts[equipment_percent < 5].sum()
df_counts = df_counts[equipment_percent >= 5]

# Ensure 'Other' category always appears in the pie chart
if 'Other' not in df_counts.columns:
    df_counts['Other'] = 0

# Add the count of 'Other' equipment types
df_counts['Other'] += other_count

df_counts = df_counts[df_counts.sum().sort_values(ascending=False).index]

# Define a color dictionary
colors=[]
#Read in color text file to get different flights to be diffrent colors on map
with open('colors.txt','r') as c_in:
	for line in c_in:
		c=str(line[0:7])
		colors.append(c)

color_dict = {equipment_type: colors[(i+3)*6 % len(colors)] for i, equipment_type in enumerate(df_counts.columns)}

# Plot the data
ax = df_counts.plot(kind='bar', stacked=True, color=[color_dict[col] for col in df_counts.columns], width=0.8)
plt.xlabel('Date')
plt.ylabel('Count')
plt.title('Occurrences of Equipment Types per Day')

# Add labels to the bars
for i in range(len(df_counts.index)):
    for j in range(len(df_counts.columns)):
        count = df_counts.iloc[i, j]
        if count > 0:
             ax.text(i, df_counts.iloc[i, :j+1].sum() - count/2, df_counts.columns[j], ha='center', va='center')

# Remove the legend
plt.legend().remove()
plt.show()


'''
import pandas as pd
import matplotlib.pyplot as plt
import glob

# Get a list of all the csv files
files = glob.glob('*_flights.csv')

# Initialize an empty list to store DataFrames
dfs = []

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
other_count = equipment_counts[equipment_percent < 0.5].sum()
equipment_counts = equipment_counts[equipment_percent >= 0.5]

# Ensure 'Other' category always appears in the pie chart
if 'Other' not in equipment_counts.index:
    equipment_counts = pd.concat([equipment_counts, pd.Series([0], index=['Other'])])

# Add the count of 'Other' equipment types
equipment_counts['Other'] += other_count

# Plot the data as a bar graph
equipment_counts.plot(kind='bar')
plt.xlabel('Equipment Type')
plt.ylabel('Count')
plt.title('Occurrences of Equipment Types')
plt.show()

# Plot the data as a pie chart
equipment_counts.plot(kind='pie', autopct='%1.1f%%')
plt.title('Occurrences of Equipment Types')
plt.show()
'''
