# Import the libraries
import numpy as np
from math import radians, sin, cos, sqrt, asin
import math
from obspy.geodetics import gps2dist_azimuth

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

print("Station 1230 distance(km) to explosion")
print("B1:", distance(63.982154, -149.121513, 63.98, -148.66))
dist = gps2dist_azimuth(63.982154, -149.121513, 63.98, -148.66)
dist_km = dist[0]/1000
print("B1:", dist_km)

print("B2:", distance(63.982154, -149.121513, 63.97, -148.68))
dist = gps2dist_azimuth(63.982154, -149.121513, 63.97, -148.68)
dist_km = dist[0]/1000
print("B2:", dist_km)

print("B3:", distance(63.4056, -148.8602, 64.01, -148.76))
dist = gps2dist_azimuth(63.4056, -148.8602, 64.01, -148.76)
dist_km = dist[0]/1000
print("B3:", dist_km) 



