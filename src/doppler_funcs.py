import os
import obspy
import json
import numpy as np
import numpy.linalg as la
import pandas as pd
import math
from pathlib import Path
from pyproj import Proj
from datetime import datetime, timezone

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

###############################################################

def make_base_dir(base_dir):

	"""
	Create a directory and its parent directories if they don't exist.

	Args:
		base_dir (str): The path of the directory to be created.

	Returns:
		None
	"""

	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path / parent
			if not current_path.exists():
				current_path.mkdir()

#######################################################################################################################

def add_wind_vector(zonal_winds, meridional_winds):
	"""
	Adds two vectors given their magnitudes and directions (angles in degrees, clockwise from the positive y-axis).

	Args:
		zonal_winds (float): Magnitude of winds eastward.
		meridional_winds (float): Magnitude of winds northward.

	Returns:
		tuple: A tuple containing the magnitude and direction (angle in degrees, clockwise from North).
	"""

	# Calculate magnitude of the resultant vector
	resultant_magnitude = math.sqrt(zonal_winds **2 + meridional_winds **2)

	# Calculate direction (angle) of the resultant vector
	resultant_angle_rad = math.atan2(meridional_winds,zonal_winds)                                 							
	resultant_angle_deg = (90 - math.degrees(resultant_angle_rad)) % 360

	return resultant_magnitude, resultant_angle_deg

################################################################################################################################

def dist_less(flight_utm_x_km, flight_utm_y_km, seismo_utm_x_km, seismo_utm_y_km):
	"""
	Checks if any seismometer is within 2 km of the flight path.

	Args:
		flight_utm_x_km (list): List of UTM x-coordinates of the flight path in kilometers.
		flight_utm_y_km (list): List of UTM y-coordinates of the flight path in kilometers.
		seismo_utm_x_km (list): List of UTM x-coordinates of the seismometers in kilometers.
		seismo_utm_y_km (list): List of UTM y-coordinates of the seismometers in kilometers.
	
	Returns:
		bool: True if any seismometer is within 2 km of the flight path, False otherwise.
	"""

	f = False
	for s in range(len(flight_utm_x_km)):
		for l in range(len(seismo_utm_x_km)):
			dist_km = np.sqrt((seismo_utm_y_km[l]-flight_utm_y_km[s])**2 +(seismo_utm_x_km[l]-flight_utm_x_km[s])**2)
			if dist_km <= 2:
				f = True
				break
			else:
				continue
	return f

#######################################################################################################################

def get_speed_of_sound(alt, closest_time, UTM_x_m, UTM_y_m):
	"""
	Calculate the speed of sound at a given altitude using atmospheric data. 
	If you want the average between a station and aircraft this function needs to be updated.

	Args:
		alt (float): Altitude in meters.
		closest_time (float): The time of the closest point on the flight path.
		UTM_x_m (float): UTM x-coordinate in meters.
		UTM_y_m (float): UTM y-coordinate in meters.

	Returns:
		float: The speed of sound in meters per second.
		float: The temperature in degrees Celsius at the given altitude.
	"""

	# Convert UTM coordinates to latitude and longitude
	lon, lat = utm_proj(UTM_x_m, UTM_y_m, inverse=True)

	input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'

	if Path(input_files).exists():
		with open(input_files, 'r') as file:
			data = json.load(file)

		# Extract data
		data_list = data['data']

		# Find the "Z" parameter and extract the value at index
		z_index = None
		hold = np.inf
		for item in data_list:
			if item['parameter'] == 'Z':
				for i in range(len(item['values'])):
					if abs(float(item['values'][i]) - float(alt/1000)) < hold:
						hold = abs(float(item['values'][i]) - float(alt/1000))
						z_index = i

		Tc = None
		for item in data_list:
			if item['parameter'] == 'T' and z_index is not None:
				Tc = -273.15 + float(item['values'][z_index])

		if Tc is not None:
			c = speed_of_sound(Tc)
		else:
			c = 311  # Default speed of sound in m/s if no temperature data is available
			Tc = -33.7
	else:
		c = 311  # Default speed of sound in m/s if no data is available
		Tc = -33.7

	return c, Tc

#######################################################################################################################

def time_check(closest_time, start_time, end_time, s_index):
	"""
	Check if the closest time is outside the start and end operating times for a given seismometer.

	Args:
		closest_time (float): The time of the closest point on the flight path.
		start_time (list): List of start times for each seismometer.
		end_time (list): List of end times for each seismometer.
		s_index (int): Index of the seismometer being checked.

	Returns:
		bool: True if the closest time is outside the start and end operating times, False otherwise.
	"""

	v = False
	# Convert the times to datetime objects
	time_airplane = datetime.fromtimestamp(closest_time, tz=timezone.utc)
	start_time_obj = datetime.strptime(start_time[s_index], "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)
	end_time_obj = datetime.strptime(end_time[s_index], "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)

	if time_airplane < start_time_obj or time_airplane > end_time_obj:
		v = True
	return v

#######################################################################################################################

def get_sta_elevation(sta):
	"""
	Get the elevation of a seismic station.

	Args:
		sta (str): The station code.

	Returns:
		float: The elevation of the station in meters.
	"""

	elev = 0
	seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/parkshwy_nodes.txt', sep="|")
	stations = seismo_data['Station']
	elevations = seismo_data['Elevation']

	for i in range(len(stations)):
		if stations[i] == sta:
			elev = elevations[i]
			break

	return elev

####################################################################################################################

def flight_lat_lon_to_utm(flight_latitudes, flight_longitudes):
	"""
	Convert flight latitude and longitude to UTM coordinates in kilometers.

	Args:
		flight_latitudes (list): List of flight latitudes in degrees.
		flight_longitudes (list): List of flight longitudes in degrees.

	Returns:
		tuple: A tuple containing lists of UTM x-coordinates and y-coordinates in kilometers,
		      and the flight path as a list of tuples containing the x and y UTM coordinates in kilometers.
	"""

	utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

	# Convert flight latitude and longitude to UTM coordinates
	flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
	flight_utm_x, flight_utm_y = zip(*flight_utm)

	# Convert UTM coordinates to kilometers
	flight_utm_x_km = [x / 1000 for x in flight_utm_x]
	flight_utm_y_km = [y / 1000 for y in flight_utm_y]
	flight_path = [(x, y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]
		
	return flight_utm_x_km, flight_utm_y_km, flight_path

#################################################################################################################################
	
def closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y):
	"""
	Calculate the closest point on a segment to a seismic station.

	Args:

		flight_utm_x1 (float): UTM x-coordinate of the first point of the flight segment.	
		flight_utm_y1 (float): UTM y-coordinate of the first point of the flight segment.
		flight_utm_x2 (float): UTM x-coordinate of the second point of the flight segment.
		flight_utm_y2 (float): UTM y-coordinate of the second point of the flight segment.
		seismo_utm_x (float): UTM x-coordinate of the seismic station.
		seismo_utm_y (float): UTM y-coordinate of the seismic station.

	Returns:
		tuple: A tuple containing the closest point on the segment, and the distance between the segment and the station.
	"""

	closest_point = None
	dist_lim = np.inf

	x = [flight_utm_x1, flight_utm_x2]
	y = [flight_utm_y1, flight_utm_y2]

	if (x[1] - x[0]) == 0:
		if (y[1]-y[0]) <= 0:
			ggg = -0.001
		else:
			ggg = 0.001
		for point in np.arange(y[0], y[1], ggg):
			xx = x[0]
			yy = point
			dist_km = np.sqrt((seismo_utm_y-yy)**2 +(seismo_utm_x-xx)**2)
			
			if dist_km < dist_lim:
				dist_lim = dist_km
				closest_point = (xx,yy)
			else:
				continue

	else: 
		m = (y[1]-y[0])/(x[1]-x[0])
		b = y[0] - m*x[0]

		if (x[1] - x[0]) <= 0:
			ggg = -0.001
		else:
			ggg = 0.001
		for point in np.arange(x[0], x[1], ggg):
			xx = point
		
			yy = m*xx + b
			dist_km = np.sqrt((seismo_utm_y-yy)**2 +(seismo_utm_x-xx)**2)
			
			if dist_km < dist_lim:
				dist_lim = dist_km
				closest_point = (xx,yy)
			else:
				continue

	return closest_point, dist_lim

#################################################################################################################################

def find_closest_point(flight_utm, seismo_utm):
	"""
	Find the closest point on a flight path to a seismic station.

	Args:
		flight_utm (list): List of UTM coordinates of the flight path.
		seismo_utm (tuple): UTM coordinates of the seismic station.

	Returns:
		tuple: A tuple containing the closest point on the flight path, 
		the distance between the flight path and the station, and the index of the closest point.
	"""

	min_distance = np.inf
	closest_point = None
	index = None
	for i in range(len(flight_utm) - 1):
		flight_utm_x1, flight_utm_y1 = flight_utm[i]
		flight_utm_x2, flight_utm_y2 = flight_utm[i + 1]
		seismo_utm_x, seismo_utm_y = seismo_utm
		point, d = closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y)
		
		if point == None:
			continue
		elif d < min_distance:
			min_distance = d
			closest_point = point
			index = i
		else:
			continue

	return closest_point, min_distance, index	

###########################################################################################################################################

def closest_time_calc(closest_p, flight_path, timestamp, index):
	"""
	Calculate the time of the closest point on the flight path to a seismic station.

	Args:
		closest_p (tuple): The closest point on the flight path in UTM coordinates, (x,y).
		flight_path (list): List of tuples containing the UTM coordinates in kilometers of the flight path, (x,y).
		timestamp (list): List of timestamps corresponding to the flight path.
		index (int): The index of the closest point, provided by flightradar24, in the flight path.

	Returns:
		float: The calculated time of the closest point.
	"""

	closest_x, closest_y = closest_p

	flight_utm_x1, flight_utm_y1 = flight_path[index]
	flight_utm_x2, flight_utm_y2 = flight_path[index + 1]

	x_timestamp_dif_vec = flight_utm_x2 - flight_utm_x1
	y_timestamp_dif_vec = flight_utm_y2 - flight_utm_y1

	cx_timestamp_dif_vec = closest_x - flight_utm_x1
	cy_timestamp_dif_vec = closest_y - flight_utm_y1

	line_vector = (x_timestamp_dif_vec, y_timestamp_dif_vec)
	cline_vector = (cx_timestamp_dif_vec, cy_timestamp_dif_vec)

	line_magnitude = np.sqrt(line_vector[0] ** 2 + line_vector[1] ** 2)
	cline_magnitude = np.sqrt(cline_vector[0] ** 2 + cline_vector[1] ** 2)

	length_ratio = cline_magnitude / line_magnitude
	closest_time = timestamp[index] + length_ratio * (timestamp[index + 1] - timestamp[index])

	return closest_time

############################################################################################################################

def avg_return(alt, speed, head, index):
	"""
	Calculate the average altitude, speed, and heading at a given index in the flight data.

	Args:
		alt (list): List of altitudes in feet.
		speed (list): List of speeds in knots.
		head (list): List of headings in degrees.
		index (int): The index in the lists for which to calculate the averages.

	Returns:
		tuple: A tuple containing the average altitude in meters, average speed in meters per second, and average heading.
	"""

	alt_avg = (alt[index] + alt[index + 1]) / 2
	alt_avg_m = alt_avg * 0.3048  # convert from feet to meters
	head_avg = (head[index] + head[index + 1]) / 2
	speed_avg = (speed[index] + speed[index + 1]) / 2
	speed_avg_mps = speed_avg * 0.514444  # convert from knots to meters/sec
	return alt_avg_m, speed_avg_mps, head_avg

#################################################################################################################

def speed_of_sound(Tc):
	"""
	Calculate the speed of sound for a given temperature.

	Parameters:
	Tc (float): Temperature in degrees celsius

	Returns:
	float: Speed of sound for temperature Tc (in m/s).
	"""

	c = 331.3+0.6*Tc
	return c

####################################################################################################################

def calc_time(t0, dist, alt, c):
	"""
	Calculate the time at which the acoustic wave reaches the station at closest approach.

	Parameters:
	t0 (float): Epoch time at which the wave is generated by the aircraft (in seconds).
	dist (float): Horizontal distance between the station and the aircraft at t0 (in meters).
	alt (float): Altitude of the aircraft at t0 (in meters) minus the elevation of the station.

	Returns:
	float: Time at which the acoustic wave reaches the station (in seconds).
	"""

	t0prime = t0 + (np.sqrt(dist**2 + alt**2))/c
	return t0prime

############################################################################################################################

def calc_ft(times, t0, f0, v0, l, c):
	"""
	Calculate the frequency at each given time using the model parameters.

	Args:
		times (list): List of time values.
		t0 (float): The time at which the central frequency of the overtones occur, 
		                when the aircraft is at the closest approach to the station.
		f0 (float): Fundamental frequency produced by the aircraft.
		v0 (float): Velocity of the aircraft.
		l (float): Distance between the station and the aircraft at the closest approach.
		c (float): Speed of sound.

	Returns:
		list: List of calculated frequency values.
	"""
	ft = []
	for tprime in times:
		t = ((tprime - t0)- np.sqrt((tprime - t0)**2-(1-v0**2/c**2)*((tprime - t0)**2-l**2/c**2)))/(1-v0**2/c**2)
		ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
								
		ft.append(ft0p)
	return ft

###################################################################################################################################################################

def S(dnew, dobs, ndata, m, mprior, cprior, tsigma):
	"""
	Calculate the data misfit using the predictions and observations.
	MISFIT FUNCTION: least squares, Tarantola (2005), Eq. 6.251
	From: https://github.com/uafgeoteach/GEOS626_seis/blob/main/hw_genlsq.ipynb

	Args:
		dnew (array): Array of predicted data.
		dobs (array): Array of observed data.
		ndata (int): Number of data points.
		m (array): Posterior model parameters.
		mprior (array): Prior model parameters.
		cprior (array): Covariance matrix for prior model.
		tsigma (float): Uncertainty in f0 measurements, Hz

	Returns:
		float: Data misfit value.
	"""

	sigma_obs = tsigma * np.ones((ndata))  # standard deviations
	cobs0 = np.diag(np.square(sigma_obs))  # diagonal covariance matrix
	m = np.array(m)
	mprior = np.array(mprior)
	Cdfac = ndata
	dnew = np.array(dnew)
	dobs = np.array(dobs)
	cobs = Cdfac * cobs0              # with normalization factor
	icobs = la.inv(cobs)
	icprior = la.inv(cprior)  # inverse covariance matrix for prior model with normalization factor

	#data misfit
	Sd = 0.5 * (dnew - dobs).T @ icobs @ (dnew - dobs)

	# model misfit (related to regularization)
	Sm = 0.5 * (m - mprior).T @ icprior @ (m - mprior)

	# total misfit
	S = Sd + Sm

	print("Model Misfit:", Sm)
	print("Data Misfit:", Sd)
	print("Total Misfit:", S)
	return Sd

###################################################################################################################################################################

def calc_f0(tprime, t0, ft0p, v0, l, c):
	"""
	Calculate the fundamental frequency produced by an aircraft (where the wave is generated) given the model parameters.

	Parameters:
	tprime (float): Time at which a frequency (ft0p) is observed on the station.
	t0 (float):  The time at which the central frequency of the overtones occur.
	ft0p (float): Frequency recorded on the seismometer, picked from the overtone doppler curve.
	v0 (float): Velocity of the aircraft.
	l (float): Distance between the station and the aircraft at the closest approach.
	c (float): Speed of sound.

	Returns:
	f0 (float): Fundamental frequency produced by the aircraft. (Frequency at the source.) 
	"""
	t = ((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2)
	f0 = ft0p*(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

	return f0

####################################################################################################################################################################################################################################################################################################################
def df(f0, v0, l, tp0, tp, c):   
    """
	Calculate the derivatives of f with respect to f0, v0, l, t0 and c.

	Parameters:
	f0 (float): Fundamental frequency produced by the aircraft.
	v0 (float): Velocity of the aircraft.
	l (float): Distance of closest approach between the station and the aircraft.
	tp0 (float): Time of that the central frequency of the overtones occur, when the aircraft is at the closest approach to the station.
	c (float): Speed of sound.
	tp (numpy.ndarray): Array of times.

	Returns:
	tuple: A tuple containing the derivatives of f with respect to f0, v0, l, t0 and c.
	"""

    #derivative with respect to f0
    f_derivef0 = (1 / (1 - (c * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))) /((c**2 - v0**2) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2))))


    #derivative of f with respect to v0
    f_derivev0 = (-f0 * v0 * (-2 * l**4 * v0**4 + l**2 * (tp - tp0)**2 * v0**6 + c**6 * (tp - tp0) * (2 * l**2 + (tp - tp0)**2 * v0**2) * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) + 
    c**2 * (4 * l**4 * v0**2 - (tp - tp0)**4 * v0**6 + l**2 * (tp - tp0) * v0**4 * (5 * tp - 5 * tp0 - 3 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) - c**4 * 
    (2 * l**4 - 3 * (tp - tp0)**3 * v0**4 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4)) - l**2 * (tp - tp0) * v0**2 * (-6 * tp + 6 * tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4)))) / (c * (c - v0) * (c + v0) * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) - c**2 * np.sqrt(l**2 + (c**4 * v0**2 * 
    (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2))
    


    #derivative of f with respect to l
    f_derivel = ((f0 * l * (tp - tp0) * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))) / 
    (c * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + 
    (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4) - 
    c**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + 
    (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2))**2))


    #derivative of f with respect to tp0
    f_derivet0 = ((f0 * l**2 * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) / 
    (c * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) - 
    c**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + 
    (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2))

    #derivative of f with respect to c
    f_derivec = (f0*v0**2*(-2*l**4*v0**4 + 2*l**2*(tp-tp0)**2*v0**6 + c**6*(tp-tp0)*(l**2 + (tp-tp0)**2*v0**2)*np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2))/c**4) + 
    c**2*(4*l**4*v0**2 - (tp-tp0)**4*v0**6 + l**2*(tp-tp0)*v0**4*(3*tp - 3*tp0 - 4*np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2))/c**4))) - 
    c**4*(l**2 + (tp-tp0)**2*v0**2)*(2*l**2 - 3*(tp-tp0)*v0**2*(-tp + tp0 + np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2))/c**4))))) / (c**2*(c - v0)*(c + v0)*
    np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2))/c**4)*np.sqrt(l**2 + (c**4*v0**2*(-tp + tp0 + np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2)) / 
    c**4))**2)/(c**2 - v0**2)**2)*(c*(-tp + tp0)*v0**2 + c*v0**2*np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2))/c**4) - 
    c**2*np.sqrt(l**2 + (c**4*v0**2*(-tp + tp0 + np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2))/c**4))**2)/(c**2 - v0**2)**2) + 
    v0**2*np.sqrt(l**2 + (c**4*v0**2*(-tp + tp0 + np.sqrt((-l**2*v0**2 + c**2*(l**2 + (tp-tp0)**2*v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2) 
    
    return f_derivef0, f_derivev0, f_derivel, f_derivet0, f_derivec


#####################################################################################################################################################################################################################################################################################################################
def invert_f(mprior, prior_sigma, coords_array, num_iterations, sigma = 10, off_diagonal = False):
	"""
	Inverts the function f using the given initial parameters and data array.

	Args:
		mprior (numpy.ndarray): Initial parameters for the function, f[0] = f0, f[1] = v0, f[2] = l, f[3] = t0, f[4] = c.
		prior_sigma (list): List of standard deviations for the prior parameters prior_sigma[0] = f0_sigma, prior_sigma[1] = v0_sigma, prior_sigma[2] = l_sigma, prior_sigma[3] = t0_sigma, prior_sigma[4] = c_sigma.
		coords_array (numpy.ndarray): Data picks along overtone doppler curve.
		num_iterations (int): Number of iterations to perform.
		sigma (float): Standard deviation for the data picks, default is 10.
		off_diagonal (bool): Whether to include off-diagonal elements in the prior covariance matrix, default is False.

	Returns:
		numpy.ndarray: The inverted parameters for the function f.
		numpy.ndarray: The covariance matrix of the posterior parameters.
		numpy.ndarray: The normalized covariance matrix of the posterior parameters.
		float: The data misfit value.
	"""

	dw,_ = coords_array.shape
	fobs = coords_array[:,1]
	tobs = coords_array[:,0]
	n = 0
 
	cprior0 = np.zeros((5,5))
	f0_sigma = prior_sigma[0]
	v0_sigma = prior_sigma[1]
	l_sigma = prior_sigma[2]
	t0_sigma = prior_sigma[3]
	c_sigma = prior_sigma[4]

	cprior0[0][0] = f0_sigma**2
	cprior0[1][1] = v0_sigma**2
	cprior0[2][2] = l_sigma**2
	cprior0[3][3] = t0_sigma**2
	cprior0[4][4] = c_sigma**2
	if off_diagonal:
		cprior0[0][3] =  -0.4*f0_sigma*t0_sigma

		cprior0[1][2] = -0.7*v0_sigma*l_sigma
		cprior0[1][4] = 0.85*v0_sigma*c_sigma

		cprior0[2][1] = -0.7*v0_sigma*l_sigma
		cprior0[2][4] = -0.7*l_sigma*c_sigma

		cprior0[3][0] =  -0.4*f0_sigma*t0_sigma

		cprior0[4][1] = 0.85*v0_sigma*c_sigma
		cprior0[4][2] = -0.7*l_sigma*c_sigma

	cprior = cprior0 * (5)

	Cd0 = np.zeros((len(fobs), len(fobs)), int)
	np.fill_diagonal(Cd0, sigma**2)
	Cd = Cd0*(dw)
	mnew = mprior.copy() #mprior is the initial guess for the parameters, mnew is the updated guess

	while n < num_iterations:
		if np.any(np.isnan(mnew)) and n == 0:
			# Handle the case where mnew contains NaN values
			return mprior, cprior0, cprior, 'Forward Model'
		elif np.any(np.isnan(mnew)):
			mnew = m
			G = G_hold
			Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
			Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
			return mnew, Cpost0, Cpost, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
		else:
			m = mnew
		f0 = m[0]
		v0 = m[1]
		l = m[2]
		t0 = m[3]
		c = m[4]

		fpred = []
		G = np.zeros((dw,5)) 

		#partial derivative matrix of f with respect to m 
		for i in range(0,dw):
			tprime = tobs[i]

			t = ((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2)
			ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
			f_derivef0, f_derivev0, f_derivel, f_derivet0, f_derivec = df(m[0], m[1], m[2], m[3], tobs[i],m[4])
			
			G[i,0:5] = [f_derivef0, f_derivev0, f_derivel, f_derivet0, f_derivec]
			fpred.append(ft0p) 
		Gm = G
		
		# steepest ascent vector (Eq. 6.307 or 6.312)
		gamma = cprior @ Gm.T @ la.inv(Cd) @ (np.array(fpred) - fobs) + (np.array(m)  - np.array(mprior)) # steepest ascent vector
		#===================================================
		# QUASI-NEWTON ALGORITHM (Eq. 6.319, nu=1)
		# approximate curvature
		H = np.identity(len(mnew)) + cprior @ Gm.T @ la.inv(Cd) @ Gm
		dm = -la.inv(H) @ gamma
		mnew = m + dm

		# Check for unreasonable parameter values
		unreasonable = (
			mnew[0] <= 5 or mnew[0] > 375 or    # f0
			mnew[1] <= 0 or mnew[1] > 350 or     # v0
			mnew[1] >= mnew[4] or  # v0 must be less than c
			mnew[2] < 0 or mnew[2] > 1e5 or      # l
			mnew[3] < 10 or mnew[3] > 240 or      # t0
			mnew[4] < 200 or mnew[4] > 400       # c
		)
		if unreasonable and n > 0:
			mnew = m
			G = G_hold
			Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
			Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
			return mnew, Cpost0, Cpost, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)

		elif unreasonable and n == 0:
			return mprior, cprior0, cprior, 'Forward Model'
		elif np.nan in mnew:
			return mprior, cprior0, cprior, 'Forward Model'
		else:
			G_hold = G.copy()
			n += 1
		print(mnew)

	Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
	Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
	F_m = S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
	del G, G_hold, Gm, H, dm, gamma, fpred, m, n, Cd, Cd0, cprior, cprior0
	return mnew, Cpost0, Cpost, F_m

#####################################################################################################################################################################################################################################################################################################################

def full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations = 4, sigma = 3, off_diagonal = False):
	"""
	Performs inversion using all picked overtones. 

	Args:
		fobs (numpy.ndarray): Picked frequency values from individual overtone inversion picks.
		tobs (numpy.ndarray): Picked time values from individual overtone inversion picks.
		peak_assos (list): List of number of peaks associated with each overtone, for indexing the fobs and tobs arrays.
		mprior (numpy.ndarray): Initial guess for the model parameters, mprior[0] = v0, mprior[1] = l, mprior[2] = t0, mprior[3] = c, mprior[4:] = f0_array.
		num_iterations (int): Number of iterations to perform for the inversion.
		sigma (float): Standard deviation for the data picks, default is 3.
		off_diagonal (bool): Whether to include off-diagonal elements in the prior covariance matrix, default is False.

	Returns:
		numpy.ndarray: The inverted parameters for the function f. Velocity of the aircraft, distance of closest approach, time of closest approach, and the fundamental frequency produced by the aircraft.
		numpy.ndarray: The covariance matrix of the inverted parameters.
		numpy.ndarray: The array of the fundamental frequency produced by the aircraft.
	"""

	w = len(mprior[4:]) #number of overtones

	qv = 0
	cprior0 = np.zeros((len(mprior),len(mprior)))

	f0_sigma = sigma_prior[0]
	v0_sigma = sigma_prior[1]
	l_sigma = sigma_prior[2]
	t0_sigma = sigma_prior[3]
	c_sigma = sigma_prior[4]

	if off_diagonal:
		cprior0[4:][2] =  -0.4*f0_sigma*t0_sigma

		cprior0[0][1] = -0.7*v0_sigma*l_sigma
		cprior0[0][3] = 0.85*v0_sigma*c_sigma

		cprior0[1][0] = -0.7*v0_sigma*l_sigma
		cprior0[1][3] = -0.7*l_sigma*c_sigma

		cprior0[2][4:] =  -0.4*f0_sigma*t0_sigma

		cprior0[3][0] = 0.85*v0_sigma*c_sigma
		cprior0[3][1] = -0.7*l_sigma*c_sigma
    
	for row in range(len(cprior0)):
		if row == 0:
			cprior0[row][row] = v0_sigma**2
		elif row == 1:
			cprior0[row][row] = l_sigma**2
		elif row == 2:
			cprior0[row][row] = t0_sigma**2
		elif row == 3:
			cprior0[row][row] = c_sigma**2
		else:
			cprior0[row][row] = f0_sigma**2
	cprior = cprior0 * (len(mprior))

	Cd0 = np.zeros((len(fobs), len(fobs)), float)
	np.fill_diagonal(Cd0, sigma**2)

	Cd = Cd0*(len(fobs))
	mnew = np.array(mprior)

	while qv < num_iterations:
		if np.any(np.isnan(mnew)) and qv == 0:
			# Handle the case where mnew contains NaN values
			return mprior, cprior0, cprior, mprior[4:], 'Forward Model'
		elif np.any(np.isnan(mnew)):
			mnew = m
			G = G_hold
			Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
			Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
			return mnew, Cpost0, Cpost, f0_array, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
		else:
			m = mnew
		v0 = m[0]
		l = m[1]
		t0 = m[2]
		c = m[3]
		f0_array = m[4:]

		fpred = []
		G = np.zeros((0,w+4))
		cum = 0
		for p in range(w):
			new_row = np.zeros(w+4)
			f0 = f0_array[p]
			
			for j in range(cum,cum+peaks_assos[p]):
				tprime = tobs[j]
				t = ((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2)
				ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

				f_derivef0, f_derivev0, f_derivel, f_derivet0, f_derivec = df(f0,v0,l,t0, tobs[j],c)
                
				new_row[0] = f_derivev0
				new_row[1] = f_derivel
				new_row[2] = f_derivet0
				new_row[3] = f_derivec
				new_row[4+p] = f_derivef0
				G = np.vstack((G, new_row))
						
				fpred.append(ft0p)
		
			cum = cum + peaks_assos[p]

		Gm = G
		
		# steepest ascent vector (Eq. 6.307 or 6.312)
		gamma = cprior @ Gm.T @ la.inv(Cd) @ (np.array(fpred) - fobs) + (np.array(m)  - np.array(mprior)) # steepest ascent vector
		#===================================================
		# QUASI-NEWTON ALGORITHM (Eq. 6.319, nu=1)
		# approximate curvature
		H = np.identity(len(mnew)) + cprior @ Gm.T @ la.inv(Cd) @ Gm
		dm = -la.inv(H) @ gamma
		mnew = m + dm

		# Check for unreasonable parameter values
		unreasonable = (
			[mn for mn in mnew[4:] if mn <= 5 or mn > 375] or   # f0
		 	mnew[0] <= 0 or mnew[0] > 350 or     # v0
		 	mnew[0] >= mnew[3] or  # v0 must be less than c
		 	mnew[1] < 0 or mnew[1] > 1e5 or      # l
		 	mnew[2] < 10 or mnew[2] > 240 or      # t0
		 	mnew[3] < 200 or mnew[3] > 400       # c
		 )

		if unreasonable and qv > 0:
			mnew = m
			G = G_hold
			Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
			Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
			return mnew, Cpost0, Cpost, f0_array, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)

		elif unreasonable and qv == 0:
			return mprior, cprior0, cprior, mprior[4:], 'Forward Model'
		elif np.nan in mnew:
			return mprior, cprior0, cprior, mprior[4:], 'Forward Model'
		else:
			# Store the current G matrix for potential rollback
			G_hold = G.copy()
		f0_array = m[4:]
		qv += 1
		print(mnew)

	Cpost = la.inv(Gm.T@la.inv(Cd)@Gm + la.inv(cprior))
	Cpost0 = la.inv(Gm.T@la.inv(Cd0)@Gm + la.inv(cprior0))
	F_m = S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
	del G, G_hold, Gm, H, dm, gamma, fpred, m, Cd, Cd0, cprior, cprior0
	return mnew, Cpost0, Cpost, mnew[4:], F_m

########################################################################################################################################################################################

def flight_list(month1, month2, first_day, last_day):
	"""
	Load flight files based on the specified months and days.

	Args:
		month1 (int): The starting month.
		month2 (int): The ending month.
		first_day (int): The first day of the range.
		last_day (int): The last day of the range.

		for only Feb use month1 = 2 and month2 = 3
		for only March use month1 = 3 and month2 = 4
		for Fab and March use month1 = 2 and month2 = 4
		for entire deployment use month1 = 2,first_day = 11,month2 = 4, and last_day = 27

	Returns:
		tuple: A tuple containing two lists - flight_files and filenames.
			   flight_files: A list of file paths for the flight files.
			   filenames: A list of filenames for the flight files.
	"""
	flight_files = []
	filenames = []

	for month in range(month1, month2):
		if month1 == 2 and month2 == 4:
			if month == 2:
				month = '02'
				for day in range(first_day, 29):
					day = str(day)
					directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
					for filename in os.listdir(directory):
						filenames.append(filename)
						f = os.path.join(directory, filename)
						if os.path.isfile(f):
							flight_files.append(f)
			elif month == 3:
				month = '03'
				for day in range(1, last_day):
					if day < 10:
						day = '0' + str(day)
						directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
						for filename in os.listdir(directory):
							filenames.append(filename)
							f = os.path.join(directory, filename)
							if os.path.isfile(f):
								flight_files.append(f)
					else:
						day = str(day)
						directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
						for filename in os.listdir(directory):
							filenames.append(filename)
							f = os.path.join(directory, filename)
							if os.path.isfile(f):
								flight_files.append(f)
		elif month1 == 2 and month2 == 3:
			month = '02'
			for day in range(first_day, last_day):
				day = str(day)
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
				for filename in os.listdir(directory):
					filenames.append(filename)
					f = os.path.join(directory, filename)
					if os.path.isfile(f):
						flight_files.append(f)
		elif month1 == 3 and month2 == 4:
			month = '03'
			for day in range(first_day, last_day):
				if day < 10:
					day = '0' + str(day)
					directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
					for filename in os.listdir(directory):
						filenames.append(filename)
						f = os.path.join(directory, filename)
						if os.path.isfile(f):
							flight_files.append(f)
				else:
					day = str(day)
					directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
					for filename in os.listdir(directory):
						filenames.append(filename)
						f = os.path.join(directory, filename)
						if os.path.isfile(f):
							flight_files.append(f)
	return flight_files, filenames

########################################################################################################################################################################################################

def load_flight_file(flight_file, filename):
	"""
	Load flight data from a CSV file and extract relevant columns.

	Args:
		flight_file (str): Path to the flight data CSV file.
		filename (str): Name of the flight file for extracting flight number and date.

	Returns:
		tuple: A tuple containing the following elements:
			- flight_utm_x_km (numpy.ndarray): UTM x coordinates of the flight path in kilometers.
			- flight_utm_y_km (numpy.ndarray): UTM y coordinates of the flight path in kilometers.
			- flight_path (numpy.ndarray): 2D array of flight path coordinates in kilometers.
			- timestamp (pandas.Series): Timestamps of the flight data.
			- alt (pandas.Series): Altitudes of the flight data.
			- speed (pandas.Series): Speeds of the flight data.
			- head (pandas.Series): Headings of the flight data.
			- flight_num (str): Flight number extracted from the filename.
			- date (str): Date extracted from the filename.
	"""
	
	flight_data = pd.read_csv(flight_file, sep=",") 
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	timestamp = flight_data['snapshot_id']  
	alt = flight_data['altitude']
	speed = flight_data['speed']
	head = flight_data['heading']
	fname = filename	
	flight_num = fname[9:18]
	date = fname[0:8]

	flight_utm_x_km, flight_utm_y_km, flight_path = flight_lat_lon_to_utm(flight_latitudes, flight_longitudes)

	return flight_utm_x_km, flight_utm_y_km, flight_path, timestamp, alt, speed, head, flight_num, date

#########################################################################################################################################################################################################

def get_equip(date, flight_num):
	"""
	Retrieve the equipment type for a specific flight number on a given date.

	Args:
		date (str): Date in the format 'YYYYMMDD'.
		flight_num (str): Flight number to search for.

	Returns:
		str: Equipment type associated with the flight number.
	"""

	equip_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_flights.csv'
	equip_data = pd.read_csv(equip_file, sep=",")
	equip_list = equip_data['equip']
	flight_list = equip_data['flight_id']

	for i_e in range(len(equip_list)):
		if str(flight_num) == str(flight_list[i_e]):
			equip = equip_list[i_e]
			break
	return equip

#########################################################################################################################################################################################################

def load_waveform(sta, start_time, spec_window=120):
	"""
	Load waveform data for a specific station and time window.

	Args:
		sta (str): Station code.
		start_time (float): Start time in seconds since the epoch.
		spec_window (int): Time window in seconds to trim the waveform data, default is 120 seconds.

	Returns:
		tuple: A tuple containing the waveform data, sampling frequency, time data correlating to waveform, and title
	"""
	
	ht = datetime.fromtimestamp(start_time + spec_window, tz=timezone.utc)

	h = ht.hour
	mins = ht.minute
	secs = ht.second
	month = ht.month
	day = ht.day

	h_u = str(h + 1)
	if h < 23:
		day2 = str(day)
		if h < 10:
			h_u = '0' + str(h + 1)
			h = '0' + str(h)
		else:
			h_u = str(h + 1)
			h = str(h)
	else:
		h_u = '00'
		day2 = str(day + 1)
	if len(str(day)) == 1:
		day = '0' + str(day)
		day2 = day

	waveform1 = "/scratch/naalexeev/NODAL/2019-0" + str(month) + "-" + str(day) + "T" + str(h) + ":00:00.000000Z.2019-0" + str(month) + "-" + str(day2) + "T" + str(h_u) + ":00:00.000000Z." + str(sta) + ".mseed"
	waveform2 = "/scratch/irseppi/500sps/2019_0" + str(month) + "_" + str(day) + "/ZE_" + str(sta) + "_DPZ.msd"
	
	if Path(waveform1).exists():
		tr = obspy.read(waveform1)
		# Trim all traces in the Stream object
		for trace in tr:
			trace.trim(trace.stats.starttime + (mins * 60) + secs - spec_window,
					trace.stats.starttime + (mins * 60) + secs + spec_window)
		data = tr[2][:]
		fs = int(tr[2].stats.sampling_rate)
		title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'
		t_wf = tr[2].times()
		if len(data) == 0:
			data = tr[1][:]
			fs = int(tr[1].stats.sampling_rate)
			title = f'{tr[1].stats.network}.{tr[1].stats.station}.{tr[1].stats.location}.{tr[1].stats.channel} − starting {tr[1].stats["starttime"]}'
			t_wf = tr[1].times()
			if len(data) == 0:
				data = tr[0][:]
				fs = int(tr[0].stats.sampling_rate)
				title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
				t_wf = tr[0].times()
		return data, fs, t_wf, title
	elif Path(waveform2).exists():
		tr = obspy.read(waveform2)
		tr[0].trim(tr[0].stats.starttime + (float(h) * 3600) + (mins * 60) + secs - spec_window,
					tr[0].stats.starttime + (float(h) * 3600) + (mins * 60) + secs + spec_window)
		data = tr[0][:]
		fs = int(tr[0].stats.sampling_rate)
		title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
		t_wf = tr[0].times()
		return data, fs, t_wf, title
	else:
		return None, None, None, None

	

#########################################################################################################################################################################################################