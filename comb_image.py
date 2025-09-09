import pandas as pd
import os
from PIL import Image, ImageDraw, ImageFont
import glob
import numpy as np
import json
from pyproj import Proj
from pathlib import Path
from src.doppler_funcs import speed_of_sound, add_wind_vector, make_base_dir
seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station_elevations = seismo_data['Elevation']
stations = seismo_data['Station']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

file_in = open('input/node_crossings_db_UTM.txt', 'r')

for line in file_in.readlines():
	text = line.split(',')
	date = text[0]
	flight_num = text[1]
	x_m =  float(text[2])  # Replace with your UTM x-coordinate
	y_m = float(text[3])  # Replace with your UTM y-coordinate
	dist_m = float(text[4])   # Distance in meters
	closest_time = float(text[5])
	alt_m = float(text[6]) 
	speed_mps = float(text[7])  # Speed in meters per second
	heading = (90 - float(text[8])) % 360


	sta = str(text[9])
	equip = text[10]
	day = str(date[6:8])
	month = str(date[4:6])

	index = None
	for i, station in enumerate(stations):
		if str(station) == str(sta):
			index = i
			break
	sta_elv = station_elevations[index]

	# Convert UTM coordinates to latitude and longitude
	lon, lat = utm_proj(x_m, y_m, inverse=True)

	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_flights.csv', sep=",")
	flight = flight_data['flight_id']
	callsign = flight_data['callsign'] 
	aircraft_id = flight_data['aircraft_id']
	for g,f_id in enumerate(flight):
		if str(f_id) == str(flight_num):
			call = callsign[g]
			id = aircraft_id[g]
			break
		else:
			continue
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/inversion_results/' + str(equip) + '_spec_c/2019-'+month+'-'+day + '/' + str(flight_num) + '/' + str(sta) + '/'
	if os.path.exists(spec_dir):
		for image in os.listdir(spec_dir):
			if not image.endswith('.png'):
				continue
			im = os.path.join(spec_dir, image)
			split_array = np.array(image.split('_'))
			plot_time = split_array[0]
	else:
		continue
	input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'
	
	if Path(input_files).exists():
		file = open(input_files, 'r')
		data = json.load(file)

		# Extract metadata
		metadata = data['metadata']
		sourcefile = metadata['sourcefile']
		datetim = metadata['time']['datetime']
		latitude = metadata['location']['latitude']
		longitude = metadata['location']['longitude']
		parameters = metadata['parameters']

		# Extract data
		data_list = data['data']

		# Convert data to a DataFrame
		data_frame = pd.DataFrame(data_list)

		# Find the "Z" parameter and extract the value at index
		z_index = None
		hold = np.inf
		for item in data_list:
			if item['parameter'] == 'Z':
				for okay in range(len(item['values'])):
					if abs(float(item['values'][okay]) - float(alt_m/1000)) < hold:
						hold = abs(float(item['values'][okay]) - float(alt_m/1000))
						z_index = okay
		for item in data_list:
			if item['parameter'] == 'T':
				Tc = -273.15 + float(item['values'][z_index])
			if item['parameter'] == 'U':
				zonal_wind = float(item['values'][z_index])
			if item['parameter'] == 'V':
				meridional_wind = float(item['values'][z_index])

		wind, az = add_wind_vector(zonal_wind, meridional_wind)
		c = speed_of_sound(Tc)
	else:
		c = 311  # Default speed of sound in m/s if no data is available
	diff = np.inf

	flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight_num) + '.csv'
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	timestamps = flight_data['snapshot_id']

	# Convert flight latitude and longitude to UTM coordinates
	flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
	flight_utm_x, flight_utm_y = zip(*flight_utm)

	# Convert UTM coordinates to kilometers
	flight_utm_x_km = [x / 1000 for x in flight_utm_x]
	flight_utm_y_km = [y / 1000 for y in flight_utm_y]
	for t in range(len(timestamps)):
		if abs(float(closest_time) - float(timestamps[t])) < diff:
			diff = abs(float(closest_time) - float(timestamps[t]))
			if t < len(flight_utm_x_km) - 1:
				direction = np.arctan2(flight_utm_y_km[t+1] - flight_utm_y_km[t], flight_utm_x_km[t+1] - flight_utm_x_km[t])
				deg = (90 -  np.degrees(direction)) % 360
			else:
				deg = heading
	dist = np.sqrt(dist_m**2 + (alt_m-sta_elv)**2)
	temp = Tc
	sound = c
	
	mnum = "FH/VT"
	font2 = ImageFont.truetype('input/Arial.ttf', 25)

			
	text1 = 'Altitude: '+str(round((alt_m-sta_elv),2))+' m\nDistance: '+str(round(dist,2))+' m\nVelocity: '+str(round(speed_mps,2))+' m/s\n               at '+str(round(deg,2))+ '\N{DEGREE SIGN}' + '\nHeading: '+str(round(heading,2))+ '\N{DEGREE SIGN}'
	text2 = 'Temperature: '+str(round(temp,1))+'\N{DEGREE SIGN}'+'C\nWind: '+str(round(wind,2))+' m/s\n         at '+str(round(az,2))+ '\N{DEGREE SIGN}\nSound Speed:\n         '+str(round(sound,2))+' m/s'
	text3 = 'Callsign: ' +  str(call) + ' (' + str(equip) + ')'

	font2 = ImageFont.truetype('input/Arial.ttf', 25)
	# Open images
	spectrogram = Image.open(im)
	# Resize images
	google_slide_width = 1280  # Width of a Google Slide in pixels
	google_slide_height = 720  # Height of a Google Slide in pixels
	# Get the path of the image file using a wildcard
	#try:
	image_path = glob.glob('/scratch/irseppi/nodal_data/plane_info/map_all_UTM/2019'+month+day+'/'+flight_num+'/'+sta+'/map_'+flight_num+'_*.png')[0]
	map_img = Image.open(image_path)
	# Only downscale if the image is larger than the target size, otherwise keep original
	target_width = int(google_slide_width * 0.28)
	#if map_img.width > target_width:
	target_height = int(target_width * map_img.height / map_img.width)
	maps = map_img.resize((target_width, target_height), Image.LANCZOS)
		#else:
		#	maps = map_img.copy()
	#except:
	#	print('No image for: ' + image_path)
	#	continue
	try:
		spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/inversion_results/' + str(equip) + '_spectrum_c/2019'+month+day+'/'+flight_num+'/'+sta+'/'+sta+'_' + str(plot_time) + '.png')
	except:
		print('No spectrum image for: ' + '/scratch/irseppi/nodal_data/plane_info/inversion_results/' + str(equip) + '_spectrum_c/2019'+month+day+'/'+flight_num+'/'+sta+'/'+sta+'_' + str(plot_time) + '.png')
		continue

	# Resize images
	google_slide_width = 1280  # Width of a Google Slide in pixels
	google_slide_height = 720  # Height of a Google Slide in pixels

	path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(equip)+'.jpg'
	if os.path.isfile(path):
		plane_img = Image.open(path)
		# 'maps' is now set above, so this line is no longer needed
		# maps = map_img.resize((int(google_slide_width *  0.28), int(google_slide_width *0.28* map_img.height / map_img.width)))
	else:
		plane_img = Image.open('hold.png')
		
	scale = 70/1280
	plane = plane_img.resize((int(google_slide_width * 0.26), int(google_slide_height * 0.26)))
	spec = spec_img.resize((int(google_slide_width * 0.31), int(google_slide_height * 0.35)))  
	#maps = map_img.resize((int(google_slide_width *  0.28), int(google_slide_width *0.28* map_img.height / map_img.width)))
	spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height)))

	# Create blank canvas
	canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

	# Paste images onto canvas
	canvas.paste(spec, (google_slide_width - spec.width+ int(spec.width/12), google_slide_height - spec.height))
	canvas.paste(maps, (google_slide_width - int(maps.width*1.05), int(plane.height)))
	canvas.paste(plane, (google_slide_width - plane.width, 0))
	canvas.paste(spectrogram, (-40, 0))
	# Draw text from files
	draw = ImageDraw.Draw(canvas)
	font = ImageFont.truetype('input/Arial.ttf', 14) 

	# Label each image
	draw.text((15, 35), '(a)', fill='black', font=font2)
	draw.text((15, 350), '(b)', fill='black', font=font2)
	draw.text((google_slide_width - int(plane.width*1.18), 7), '(c)', fill='black', font=font2)
	draw.text((google_slide_width - int(plane.width*1.18), int(plane.height) + int(plane.height*0.05)), '(d)', fill='black', font=font2)
	draw.text((google_slide_width - int(plane.width*1.14), google_slide_height - spec.height + 20), '(e)', fill='black', font=font2)

	draw.text((google_slide_width - 305, 405), text1, fill='black', font=font)			
	draw.text((google_slide_width - 155, 405), text2,fill='black', font=font)
	bbox = draw.textbbox((google_slide_width - plane.width, 0), text3, font=font)
	draw.rectangle(bbox, fill="white")
	draw.text((google_slide_width - plane.width, 0), text3, fill='black', font=font)

	BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inverse_final_database/'
	make_base_dir(BASE_DIR)
	name= BASE_DIR +str(equip)+'_'+ '2019'+month+day+'_'+str(flight_num)+'_' + str(closest_time) + '_' + str(sta) + '_' + str(equip)+'.pdf'

	# Save combined image
	canvas.save(name, 'PDF', resolution=600.0)

file_in.close()
