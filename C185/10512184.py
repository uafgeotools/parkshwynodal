import pandas as pd
import os
from PIL import Image, ImageDraw, ImageFont
import glob
import numpy as np
import json
from pyproj import Proj
from doppler_funcs import speed_of_sound, add_vectors, make_base_dir

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

day = str(21)
month = '02'

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/529754214_list.txt','r')

flight_id = '529754214'	
aircraft_id	 = '1051218'	
equip = 'C185'	
callsign = 'N125KT'	

spec_dir = '/scratch/irseppi/nodal_data/plane_info/C185_spec_c/2019-02-21/529754214'
flight_file = '/scratch/irseppi/nodal_data/flightradar24/20190221_positions/20190221_' + flight_id + '.csv'
flight_data = pd.read_csv(flight_file, sep=",")
flight_latitudes = flight_data['latitude']
flight_longitudes = flight_data['longitude']

# Convert flight latitude and longitude to UTM coordinates
flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
flight_utm_x, flight_utm_y = zip(*flight_utm)

# Convert UTM coordinates to kilometers
flight_utm_x_km = [x / 1000 for x in flight_utm_x]
flight_utm_y_km = [y / 1000 for y in flight_utm_y]
times = flight_data['snapshot_id']

alt_list = []
lon_list = []
lat_list = []
sta_list = []
time_list = []
dist_list = []
vel_list = []
head_list = []

for line in file_in.readlines():
	#date,flight_num,closest_x_m,closest_y_m,dist_m,closest_time,alt_avg_m,speed_avg_mps,head_avg,station
	text = line.split(',')
	time_list.append(float(text[5]))
	sta_list.append(str(text[9]))
	dist_list.append(float(text[4]))
	vel_list.append(float(text[7]))
	head_list.append(float(text[8]))
	alt_list.append(float(text[6])) #convert between feet and km

	x =  float(text[2])  # Replace with your UTM x-coordinate
	y = float(text[3])  # Replace with your UTM y-coordinate

	# Convert UTM coordinates to latitude and longitude
	lon, lat = utm_proj(x, y, inverse=True)
	lat_list.append(lat)
	lon_list.append(lon)
file_in.close()

# Convert the lists to numpy arrays
lat_list = np.array(lat_list)
lon_list = np.array(lon_list)
alt_list = np.array(alt_list)
sta_list = np.array(sta_list)
time_list = np.array(time_list)
head_list = np.array(head_list)
alt_list = np.array(alt_list)
vel_list = np.array(vel_list)

for station in os.listdir(spec_dir):
	sta = os.path.join(spec_dir, station)
	for image in os.listdir(sta):
		im = os.path.join(sta, image)
		split_array = np.array(image.split('_'))
		time = str(split_array[0])
		for ggy, sta in enumerate(sta_list):
			if sta == station and ggy!= 15:
				lat = lat_list[ggy]
				lon = lon_list[ggy]
				alt = alt_list[ggy]
				t_db = time_list[ggy]
				speed = vel_list[ggy]
				head = head_list[ggy]
				dist_g = dist_list[ggy]
				id = ggy
				break

	input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(t_db) + '_' + str(lat) + '_' + str(lon) + '.dat'
	file =  open(input_files, 'r') 
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
				if abs(float(item['values'][okay]) - float(alt/1000)) < hold:
					hold = abs(float(item['values'][okay]) - float(alt/1000))
					z_index = okay
	for item in data_list:
		if item['parameter'] == 'T':
			Tc = - 273.15 + float(item['values'][z_index])
		if item['parameter'] == 'U':
			zonal_wind = float(item['values'][z_index])
		if item['parameter'] == 'V':
			meridional_wind = float(item['values'][z_index])

	if zonal_wind > 0:
		v1_angle = 90
	else:
		v1_angle = 270
	if meridional_wind > 0:
		v2_angle = 0
	else:
		v2_angle = 180
	wind, az = add_vectors(zonal_wind, v1_angle, meridional_wind, v2_angle)
	c = speed_of_sound(Tc)
	diff = np.inf
	for t in range(len(times)):
		if abs(float(time) - float(times[t])) < diff:
			diff = abs(float(time) - float(times[t]))
			direction = np.arctan2(flight_utm_y_km[t+1] - flight_utm_y_km[t], flight_utm_x_km[t+1] - flight_utm_x_km[t])
		else:
			continue

	deg = (90 -  np.degrees(direction)) % 360
	dist = np.sqrt(dist_g**2 + alt**2)
	temp = Tc
	sound = c
	
	mnum = "FH/VT"
	font2 = ImageFont.truetype('input/Arial.ttf', 25)
	diff = np.inf
			
	text1 = 'Altitude: '+str(round(alt,2))+' m\nDistance: '+str(round(dist,2))+' m\nVelocity: '+str(round(speed,2))+' m/s\n               at '+str(round(deg,2))+ '\N{DEGREE SIGN}' + '\nHeading: '+str(round(head,2))+ '\N{DEGREE SIGN}'
	text2 = 'Temperature: '+str(round(temp,1))+'\N{DEGREE SIGN}'+'C\nWind: '+str(round(wind,2))+' m/s\n         at '+str(round(az,2))+ '\N{DEGREE SIGN}\nSound Speed:\n         '+str(round(sound,2))+' m/s'
	text3 = 'Callsign: ' +  str(callsign) + ' (' + str(equip) + ')'

	font2 = ImageFont.truetype('input/Arial.ttf', 25)
	# Open images
	spectrogram = Image.open(im)

	# Get the path of the image file using a wildcard
	image_path = glob.glob('/scratch/irseppi/nodal_data/plane_info/map_all_UTM/20190221/'+flight_id+'/'+station+'/map_'+flight_id+'_*')[0]

	map_img = Image.open(image_path)
	spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/C185_spectrum_c/20190221/'+flight_id+'/'+station+'/'+station+'_' + str(time) + '.png')

	# Resize images
	google_slide_width = 1280  # Width of a Google Slide in pixels
	google_slide_height = 720  # Height of a Google Slide in pixels

	path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(equip)+'.jpg'
	if os.path.isfile(path):
		plane_img = Image.open(path)
		
	else:
		plane_img = Image.open('hold.png')
		
	scale = 70/1280
	plane = plane_img.resize((int(google_slide_width * 0.26), int(google_slide_height * 0.26)))
	spec = spec_img.resize((int(google_slide_width * 0.31), int(google_slide_height * 0.35)))  
	maps = map_img.resize((int(google_slide_width *  0.28), int(google_slide_width *0.28* map_img.height / map_img.width)))
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
	#draw.text((google_slide_width - int(plane.width*1.5), google_slide_height - spec.height - spec.height/2), '['+str(mnum)+']', fill='black', font=font2)
	draw.text((15, 350), '(b)', fill='black', font=font2)
	draw.text((google_slide_width - int(plane.width*1.18), 7), '(c)', fill='black', font=font2)
	draw.text((google_slide_width - int(plane.width*1.18), int(plane.height) + int(plane.height*0.05)), '(d)', fill='black', font=font2)
	draw.text((google_slide_width - int(plane.width*1.14), google_slide_height - spec.height + 20), '(e)', fill='black', font=font2)

	draw.text((google_slide_width - 305, 405), text1, fill='black', font=font)			
	draw.text((google_slide_width - 155, 405), text2,fill='black', font=font)
	bbox = draw.textbbox((google_slide_width - plane.width, 0), text3, font=font)
	draw.rectangle(bbox, fill="white")
	draw.text((google_slide_width - plane.width, 0), text3, fill='black', font=font)
	#show image

	BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/fig_paper/'
	make_base_dir(BASE_DIR)
	name= BASE_DIR + '20190221_'+(flight_id)+'_'+time+'_'+str(station)+'_'+str(equip)+'.png'

	# Save combined image
	canvas.save(name)
							
