# Import the libraries
import pandas as pd
import numpy as np
from math import radians, sin, cos, sqrt, asin
import geopandas as gpd 
import fiona
from shapely.geometry import MultiLineString
import shapely.wkt
from shapely.geometry import Point
from geopandas import GeoDataFrame
import matplotlib.pyplot as plt
import math
import pygmt

# Define a function to calculate the distance between two points given their latitudes and longitudes
def distance(lat1, lon1, lat2, lon2):
  # Convert degrees to radians
  lat1 = math.radians(lat1)
  lon1 = math.radians(lon1)
  lat2 = math.radians(lat2)
  lon2 = math.radians(lon2)

  # Apply the haversine formula
  dlon = lon2 - lon1
  dlat = lat2 - lat1
  a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
  c = 2 * math.asin(math.sqrt(a))
  r = 6371 # Radius of the earth in km
  return c * r # Distance in km

# Define a function to take two different size lists of lat and lon values and compare each set of lat and lon values to all the lat and lon values in the other list and then find the shortest distance between those two points and return their index numbers
def compare_lists(list1, list2):
  # Initialize the minimum distance and the index numbers
  color=[]
  index1 = -1
  index2 = -1
  
  fig = pygmt.Figure()
  fig.basemap(region=[-151,-148, 62.2, 64.6], projection="M20c", frame="a")
  fig.coast(land="white", water="skyblue")
  fig.plot(data ='Alaska_Railroad.txt', style="c0.1c", fill='black', label='Train Tracks', connection='f')
  my_labels = {'red' : '> 5km from track', 'orange' : '> 4km from track (< 5km)', 'yellow' : '> 3km from track (< 4km)', 'green' : '> 2km from track (< 3km)', 'blue' : '> 1km from track (< 2km)', 'purple' : '< 1km from track'}
  
  # Loop through the first list
  for i in range(0,302):
    # Get the latitude and longitude of the current point in the first list
    lat1 = list1[i][0]
    lon1 = list1[i][1]
    min_dist = float("inf")
    
    # Loop through the second list
    for j in range(len(list2)):
      # Get the latitude and longitude of the current point in the second list
      lat2 = list2[j][0]
      lon2 = list2[j][1]
      
      # Calculate the distance between the two points
      dist = distance(lat1, lon1, lat2, lon2)
      
      # Update the minimum distance and the index numbers if the current distance is smaller than the previous minimum distance
      if dist < min_dist and i>1:
        min_dist = dist
        index1 = i
        index2 = j
       
    #Print distaces at certain station to the given (closest) track lat/lon
    n = str(gdf1["Name"][i])
    nr = n.replace('ZE.', '')
    nr = nr.replace(' - 1', ' ')
    n = nr.split(" ")
    print(n[0],"-", round(min_dist, 4), "km", list1[i])
    
    if min_dist > 5:
      color='red'
      l='> 5km from track'
    if min_dist > 4 and min_dist < 5:
      color='orange'
      l='> 4km from track (< 5km)'
    if min_dist > 3 and min_dist < 4:
      color='yellow'
      l='> 3km from track (< 4km)'
    if min_dist > 2 and min_dist < 3:
      color='green'
      l='> 2km from track (< 3km)'
    if min_dist > 1 and min_dist < 2:
      color='blue'
      l='> 1km from track (< 2km)'
    if min_dist > 0 and min_dist < 1:
      color='purple'
      l='< 1km from track'
      
    fig.plot(x=list1[i][1], y=list1[i][0], style="c0.1c", fill=color, label=my_labels[color])
    my_labels[color] = None

  fig.legend()
  fig.show()
  
# Return the index numbers of the two points with the shortest distance
  return index1, index2, min_dist

# IMPORT KML DRIVER
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'

# READ KML file to a geopandas dataframe 
geo_df1 = gpd.read_file('ZE_NODAL.kml',driver='LIBKML')
geo_df2 = pd.read_csv('Alaska_Railroad.txt', sep=r"\s{2,}", engine="python", header=None)

# Create Pandas Dataframe from GeoPandas 
df1= pd.DataFrame(geo_df1)
df2= pd.DataFrame(geo_df2)

# Create GeoDataframe from GeoPandas
gdf1 = GeoDataFrame(df1) 
gdf2 = GeoDataFrame(df2)

gdf1 = gdf1.set_geometry("geometry")
boroughs_4326 = gdf1.to_crs("EPSG:4326")
boroughs_4326.crs

#Extract lat and lon
gdf1['lon'] = gdf1["geometry"].x
gdf1['lat'] = gdf1["geometry"].y

gdf2['lon'] = gdf2[0]
gdf2['lat'] = gdf2[1]

# Extract the coordinates as arrays
P1 = np.array([gdf1['lat'].to_numpy(),gdf1['lon'].to_numpy()])
P2 = np.array([gdf2['lat'].to_numpy(),gdf2['lon'].to_numpy()])

P1=P1.T
P2=P2.T

# Call the function and print the result
result = compare_lists(P1, P2)
print(result)
i=result[0]
print('Station:',gdf1["Name"][i])
print(f'The shortest distance is {result[2]:.3f} km between point {result[0]:.3f} in ZE_NODAL.kml and point {result[1]:.3f} in Alaska_Railroad.kml.')

