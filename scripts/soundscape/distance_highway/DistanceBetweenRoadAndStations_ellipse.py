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
from obspy.geodetics import gps2dist_azimuth

# Define a function to take two different size lists of lat and lon values and compare each set of lat and lon values to all the lat and lon values in the other list and then find the shortest distance between those two points and return their index numbers
def compare_lists(list1, list2):

  # Initialize the minimum distance and the index numbers
  color=[]
  index1 = -1
  index2 = -1
  tmin = float("inf")

  fig = pygmt.Figure()
  fig.basemap(region=[-150.5,-148.5, 62.2, 64.6], projection="M20c", frame="a")
  fig.coast(land="white", water="skyblue")
  fig.plot(data ='parks_highway.txt', style="c0.1c", fill='black', label='Road', connection='f')
  my_labels = {'red' : '> 0.05km from road', 'orange' : '> 0.04km from road (< 0.05km)', 'yellow' : '> 0.03km from road (< 0.04km)', 'green' : '> 0.02km from road (< 0.03km)', 'blue' : '> 0.01km from road (< 0.02km)', 'purple' : '< 0.01km from road'}

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
      dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
      dist_km = dist[0]/1000

      # Update the minimum distance and the index numbers if the current distance is smaller than the previous minimum distance
      if dist_km < min_dist and i>1:
        min_dist = dist_km
        if min_dist < tmin:
          tmin = min_dist
          index1 = i
          index2 = j  

    print(min_dist, 'km - ZE_NODAL index:', i, gdf1["Name"][i], list1[i])
    if min_dist > 0.05:
      color='red'
      l='> 0.05km from road'
    if min_dist > 0.04 and min_dist < 0.05:
      color='orange'
      l='> 0.04km from road (< 0.05km)'
    if min_dist > 0.03 and min_dist < 0.04:
      color='yellow'
      l='> 0.03km from road (< 0.04km)'
    if min_dist > 0.02 and min_dist < 0.03:
      color='green'
      l='> 0.02km from road (< 0.03km)'
    if min_dist > 0.01 and min_dist < 0.02:
      color='blue'
      l='> 0.01km from road (< 0.02km)'
    if min_dist > 0.00 and min_dist < 0.01:
      color='purple'
      l='< 0.01km from road'
    
    fig.plot(x=list1[i][1], y=list1[i][0], style="c0.1c", fill=color, label=my_labels[color])
    my_labels[color] = None
    
  fig.legend()
  fig.show()
  
# Return the index numbers of the two points with the shortest distance
  return index1, index2, tmin

# IMPORT KML DRIVER
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'

# READ KML file to a geopandas dataframe 
geo_df1 = gpd.read_file('ZE_NODAL.kml',driver='LIBKML')
geo_df2 = pd.read_csv('parks_highway.txt', sep=r",", engine="python", header=None)

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

print(f'The shortest distance is {result[2]:.3f} km between point {result[0]:.3f} in ZE_NODAL.kml and point {result[1]:.3f} in parks_highway.kml.')


