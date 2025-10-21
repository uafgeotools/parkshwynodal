import faulthandler

from obspy import read
from pathlib import Path

faulthandler.enable()

#Look for directory in 50sps folder and if it doesn't exist create it
def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()
#loop over the two months the data were running
for month in range(2,4):
	#for Febuary use days 11-29
	if month == 2:
		for day in range (11,29):
			#call the make base directory function
			BASE_DIR = "/home/irseppi/nodal_data/50sps/2019_0"+str(month)+"_"+str(day)
			make_base_dir(BASE_DIR)
			
			#looking for the stations 1001-1306
			for z in range(0,4):
				for y in range(0,10):
					#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
					try:
						#reading in 500sps data
						tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_1"+str(z)+str(y)+"*.msd")
						print('READ')
						
						#copying and resampling data to 50sps
						tr_new = tr.copy()
						tr_new.resample(50)
						
						#going through each file of the resampled data to create a file name and writing this to be saved in the 50sps folder
						for x in range(len(tr_new)):
							save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
							print(save_name)

							save_path = (Path(BASE_DIR) / save_name)

							tr_new[x].write(str(save_path), format='MSEED')
							print('SAVED')
					except:
						continue
			#for stations 1501-1589
			for y in range(0,9):
				#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
				try:
					#reading in 500sps data
					tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_15"+str(y)+"*.msd")
					print('READ')
					
					#copying and resampling data to 50sps
					tr_new = tr.copy()
					tr_new.resample(50)

					#going through each file of the resampled data to create a file name and writing this to be saved in the 50sps folder
					for x in range(len(tr_new)):
						save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
						print(save_name)

						save_path = (Path(BASE_DIR) / save_name)

						tr_new[x].write(str(save_path), format='MSEED')
						print('SAVED')
				except:
					continue

			#Looking for stations starting with 5
			#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
			try:
				#reading in 500sps data
				tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_5*.msd")
				print('READ')

				#copying and resampling data to 50sps
				tr_new = tr.copy()
				tr_new.resample(50)

				#going through each file of the resampled data to create a file name and writing this to be saved in the 50sps folder
				for x in range(len(tr_new)):
					save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
					print(save_name)

					save_path = (Path(BASE_DIR) / save_name)

					tr_new[x].write(str(save_path), format='MSEED')  
					print('SAVED')
			except:
				continue
				
			#looking for stations starting with 9
			#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
			try:
				#reading in 500sps data
				ts = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_9*.msd")
				print('READ')

				#copying and resampling data to 50sps
				ts_new = ts.copy()
				ts_new.resample(50)

				#going through each file of the resampled data to create a file name and writing this to be saved in the 50sps folder
				for x in range(len(ts_new)):
					save_name = BASE_DIR + "/ZE_{}_{}.msd".format(ts_new[x].stats.station,ts_new[x].stats.channel)
					print(save_name)

					save_path = (Path(BASE_DIR) / save_name)

					ts_new[x].write(str(save_path), format='MSEED')  
					print('SAVED')
			except:
				continue
	#for March use days 1-28
	else: 
		for day in range (1,10):
			BASE_DIR = "/home/irseppi/nodal_data/50sps/2019_0"+str(month)+"_0"+str(day)
			make_base_dir(BASE_DIR)
			
			#Looking for the stations 1001-1306
			for z in range(0,4):
				for y in range(0,10):
					#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
					try:
						#reading in 500sps data
						tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_0"+str(day)+"/ZE_1"+str(z)+str(y)+"*.msd")
						print('READ')

						#copying and resampling data to 50sps
						tr_new = tr.copy()
						tr_new.resample(50)
			    
						#going through each file of the resampled data to create a file name and writing this to be saved in the 50sps folder
						for x in range(len(tr_new)):
							save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
							print(save_name)

							save_path = (Path(BASE_DIR) / save_name)

							tr_new[x].write(str(save_path), format='MSEED')  
							print('SAVED')
					except:
						continue
			#for stations 1501-1589
			for y in range(0,9):
				#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
				try:
					tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_0"+str(day)+"/ZE_15"+str(y)+"*.msd")
					print('READ')

					tr_new = tr.copy()
					tr_new.resample(50)
					
					for x in range(len(tr_new)):
						save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
						print(save_name)

						save_path = (Path(BASE_DIR) / save_name)

						tr_new[x].write(str(save_path), format='MSEED')  
						print('SAVED')
				except:
					continue
			#looking for stations starting with 5
			#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
			try:
				tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_0"+str(day)+"/ZE_5*.msd")
				print('READ')

				tr_new = tr.copy()

				tr_new.resample(50)

				for x in range(len(tr_new)):
					save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
					print(save_name)

					save_path = (Path(BASE_DIR) / save_name)

					tr_new[x].write(str(save_path), format='MSEED')  
					print('SAVED')

			except:
				continue
			#looking for stations starting with 9
			#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
			try:
				ts = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_0"+str(day)+"/ZE_9*.msd")
				print('READ')

				ts_new = ts.copy()

				ts_new.resample(50)

				for x in range(len(ts_new)):
					save_name = BASE_DIR + "/ZE_{}_{}.msd".format(ts_new[x].stats.station,ts_new[x].stats.channel)
					print(save_name)

					save_path = (Path(BASE_DIR) / save_name)

					ts_new[x].write(str(save_path), format='MSEED')  
					print('SAVED')

			except:
				continue


		for day in range (10,27):
			BASE_DIR = "/home/irseppi/nodal_data/50sps/2019_0"+str(month)+"_"+str(day)
			make_base_dir(BASE_DIR)
			
			#for the stations 1001-1306
			for z in range(0,4):
				for y in range(0,10):
					#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
					try:
						tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_1"+str(z)+str(y)+"*.msd")
						print('READ')

						tr_new = tr.copy()
						tr_new.resample(50)
			    
						for x in range(len(tr_new)):
							save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
							print(save_name)

							save_path = (Path(BASE_DIR) / save_name)

							tr_new[x].write(str(save_path), format='MSEED')  
							print('SAVED')
					except:
						continue
			#for stations 1501-1589
			for y in range(0,9):
				#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
				try:
					tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_15"+str(y)+"*.msd")
					print('READ')

					tr_new = tr.copy()
					tr_new.resample(50)
		    
					for x in range(len(tr_new)):
						save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
						print(save_name)

						save_path = (Path(BASE_DIR) / save_name)

						tr_new[x].write(str(save_path), format='MSEED')  
						print('SAVED')
				except:
					continue
			#looking for stations starting with 5
			#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
			try:
				tr = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_5*.msd")
				print('READ')

				tr_new = tr.copy()

				tr_new.resample(50)

				for x in range(len(tr_new)):
					save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
					print(save_name)

					save_path = (Path(BASE_DIR) / save_name)

					tr_new[x].write(str(save_path), format='MSEED')  
					print('SAVED')
			except:
				continue
			
			#looking for stations starting with 9
			#see if there is a 500sps file path with date and station. If it doesn't exist continue through loop
			try:
				ts = read("/home/irseppi/nodal_data/500sps/2019_0"+str(month)+"_"+str(day)+"/ZE_9*.msd")
				print('READ')

				ts_new = ts.copy()

				ts_new.resample(50)

				for x in range(len(tr_new)):
					save_name = BASE_DIR + "/ZE_{}_{}.msd".format(ts_new[x].stats.station,ts_new[x].stats.channel)
					print(save_name)

					save_path = (Path(BASE_DIR) / save_name)

					ts_new[x].write(str(save_path), format='MSEED')  
					print('SAVED')
			except:
				continue
					
