import pandas as pd
import os
from PIL import Image, ImageDraw, ImageFont
import glob
import numpy as np
import json
from pyproj import Proj, Geod
from pathlib import Path
from src.doppler_funcs import speed_of_sound, add_wind_vector, make_base_dir
from pdf2image import convert_from_path

def load_pdf_as_image(pdf_path, page=0, dpi=600):
	# Convert first page of PDF to PIL Image if it's a PDF, otherwise open as image
	if str(pdf_path).lower().endswith('.pdf'):
		images = convert_from_path(pdf_path, dpi=dpi)
		return images[page]
	else:
		print("Not a PDF file for conversion:", pdf_path)

seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station_elevations = seismo_data['Elevation']
stations = seismo_data['Station']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')
geod = Geod(ellps='WGS84')
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt', 'r')

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
	file_check = '/scratch/irseppi/nodal_data/plane_info/inverse_final_database_NGT/' + str(equip)+'_'+ '2019'+month+day+'_'+str(flight_num)+'_' + str(closest_time) + '_' + str(sta) + '_' + str(equip)+'.pdf'
	if Path(file_check).exists():
		continue

	index = None
	for i, station in enumerate(stations):
		if str(station) == str(sta):
			index = i
			break
	sta_elv = station_elevations[index]
	sta_lat = seismo_latitudes[index]
	sta_lon = seismo_longitudes[index]
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
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt/' + str(equip) + '_spec_c/2019-'+month+'-'+day + '/' + str(flight_num) + '/' + str(sta) + '/'
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
	
	try:
		file = open(input_files, 'r')
		data = json.load(file)

		# Extract metadata
		metadata = data['metadata']
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
				Tc_air = -273.15 + float(item['values'][z_index])
			if item['parameter'] == 'U':
				zonal_wind_air = float(item['values'][z_index])
			if item['parameter'] == 'V':
				meridional_wind_air = float(item['values'][z_index])
		c_air = speed_of_sound(Tc_air)
		file.close()
		input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data_nodes/' + str(closest_time) + '_' + str(sta_lat) + '_' + str(sta_lon) + '.dat'
		file = open(input_files, 'r')
		data = json.load(file)

		# Extract metadata
		metadata = data['metadata']
		parameters = metadata['parameters']

		# Extract data
		data_list = data['data']

		# Convert data to a DataFrame
		data_frame = pd.DataFrame(data_list)

		for item in data_list:
			# Find the "Z" parameter and extract the value at index
			z_index = None
			hold = np.inf
			for item in data_list:
				if item['parameter'] == 'Z0':
					ground_height = float(item['values'][0])
					break 

			for item in data_list:
				if item['parameter'] == 'Z':
					for i in range(len(item['values'])):
						if abs(float(item['values'][i]) - float(ground_height)) < hold:
							hold = abs(float(item['values'][i]) - float(ground_height))
							z_index = i
			for item in data_list:
				if item['parameter'] == 'T':
					Tc_sta = - 273.15 + float(item['values'][z_index])
				if item['parameter'] == 'V':
					meridional_wind_sta = float(item['values'][z_index])
				if item['parameter'] == 'U':
					zonal_wind_sta = float(item['values'][z_index])
		wind_sta, az_sta = add_wind_vector(zonal_wind_sta, meridional_wind_sta)
		c_sta = speed_of_sound(Tc_sta)
		file.close()
		c = (c_air + c_sta) / 2
		Tc = (Tc_air + Tc_sta) / 2
		zonal_wind = zonal_wind_air #(zonal_wind_air + zonal_wind_sta) / 2
		meridional_wind = meridional_wind_air #(meridional_wind_air + meridional_wind_sta) / 2
		wind, az = add_wind_vector(zonal_wind, meridional_wind)
	except:
		c = 320  # Default speed of sound in m/s if no data is available
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
	dist = round(dist,2)
	alt = round((alt_m-sta_elv),2)
	if alt_m == 0:
		dist = "--"
		alt = "--"
	temp = Tc
	sound = c
	

	_, backazimuth, _ = geod.inv(lon, lat, seismo_longitudes[index], seismo_latitudes[index])

	text1 = 'Altitude: '+str(alt)+' m\nDistance: '+str(dist)+' m\n               at '+str(round(backazimuth,2))+ '\N{DEGREE SIGN}\nVelocity: '+str(round(speed_mps,2))+' m/s\n               at '+str(round(deg,2))+ '\N{DEGREE SIGN}'
	text2 = 'Temperature: '+str(round(temp,1))+'\N{DEGREE SIGN}'+'C\nWind: '+str(round(wind,2))+' m/s\n         at '+str(round(az,2))+ '\N{DEGREE SIGN}\nSound Speed:\n         '+str(round(sound,2))+' m/s'
	text3 = 'Callsign: ' +  str(call) + ' (' + str(equip) + ')'

	font2 = ImageFont.truetype('input/fig_style/Arial.ttf', (25/96)*600)  # Adjust size for 600 DPI


	# Get the path of the image file using a wildcard
	image_path = glob.glob('/scratch/irseppi/nodal_data/plane_info/map_all_UTM/2019'+month+day+'/'+flight_num+'/'+sta+'/map_'+flight_num+'_*.pdf')[0]
 
		
	spectrogram = Image.open(im)

	map_img = load_pdf_as_image(image_path)
	
		
	spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/inversion_results_ngt/' + str(equip) + '_spectrum_c/2019'+month+day+'/'+flight_num+'/'+sta+'/'+sta+'_' + str(plot_time) + '.png')

	path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(equip)+'.jpg'
	if os.path.isfile(path):
		plane_img = Image.open(path)
	else:
		plane_img = Image.open('/home/irseppi/REPOSITORIES/parkshwynodal/input/hold.png')
		
	# Resize images
	google_slide_width = 1280  # Width of a Google Slide in pixels
	google_slide_height = 720  # Height of a Google Slide in pixels

	# For example, for an 8x5 inch PDF at 600 DPI:
	canvas_width = 8000  # 1280/96 * 600
	canvas_height = 4500 # 720/96 * 600

	canvas = Image.new('RGB', (canvas_width, canvas_height), 'white')

	# Resize images to fit the new canvas, keeping their quality
	plane = plane_img.resize((int(canvas_width * 0.26), int(canvas_height * 0.26)), Image.LANCZOS)
	spec = spec_img.resize((int(canvas_width * 0.31), int(canvas_height * 0.35)), Image.LANCZOS)
	maps = map_img.resize((int(canvas_width *  0.28), int(canvas_width *0.28* map_img.height / map_img.width)), Image.LANCZOS)
	spectrogram = spectrogram.resize((int(canvas_width * 0.73), int(canvas_height)), Image.LANCZOS)


	# Paste images onto canvas (positions and sizes adjusted for 600 DPI)
	canvas.paste(spectrogram, (-int(20 / 96 * 600), 0))  # -20 px at 96 DPI scaled to 600 DPI
	canvas.paste(spec, (canvas_width - spec.width + int(spec.width / 12), canvas_height - spec.height))
	canvas.paste(maps, (canvas_width - int(maps.width * 1.05), int(plane.height)))
	canvas.paste(plane, (canvas_width - plane.width, 0))
	# Draw text from files
	draw = ImageDraw.Draw(canvas)
	font = ImageFont.truetype('input/fig_style/Arial.ttf', (15/96)*600)  # Adjust size for 600 DPI

	# Label each image (adjust positions and font sizes for 600 DPI)
	label_font_size = int((25/96)*600)
	label_font = ImageFont.truetype('input/fig_style/Arial.ttf', label_font_size)

	# Example y-offsets for labels, scaled for DPI
	draw.text((int(15/96*600), int(35/96*600)), '(a)', fill='black', font=label_font)
	draw.text((int(15/96*600), int(350/96*600)), '(b)', fill='black', font=label_font)
	draw.text((canvas_width - int(plane.width*1.18), int(7/96*600)), '(c)', fill='black', font=label_font)
	draw.text((canvas_width - int(plane.width*1.18), int(plane.height) + int(plane.height*0.05)), '(d)', fill='black', font=label_font)
	draw.text((canvas_width - int(plane.width*1.14), canvas_height - spec.height + int(20/96*600)), '(e)', fill='black', font=label_font)

	# Adjust text box positions for DPI
	draw.text((canvas_width - int(305/96*600), int(412/96*600)), text1, fill='black', font=font)			
	draw.text((canvas_width - int(155/96*600), int(412/96*600)), text2, fill='black', font=font)
	bbox = draw.textbbox((canvas_width - plane.width, 0), text3, font=font)
	draw.rectangle(bbox, fill="white")
	draw.text((canvas_width - plane.width, 0), text3, fill='black', font=font)

	BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inverse_final_database_NGT/'
	make_base_dir(BASE_DIR)
	name= BASE_DIR +str(equip)+'_'+ '2019'+month+day+'_'+str(flight_num)+'_' + str(closest_time) + '_' + str(sta) + '_' + str(equip)+'.pdf'

	# Save as PDF
	canvas.save(name, "PDF", resolution=700.0)
file_in.close()
