import numpy as np
from scipy.signal import spectrogram
import obspy
import numpy as np
import os
import obspy
import datetime
from pathlib import Path
from scipy import signal

text = open('input/all_station_crossing_db.txt', 'r')
#output = open('okay.txt','w')

for line in text.readlines():
	val = line.split(',')
	time = val[2]
	ht = datetime.datetime.utcfromtimestamp(int(val[2]))
	mins = ht.minute
	secs = ht.second
	h = ht.hour
	month1 = ht.month
	day1 = ht.day
	station = str(val[5])
	plane = str(val[6])
	tim = 120
	
	if ht.month == 2 and ht.day >= 11 or ht.month == 3:
		day_of_year = str((ht - datetime.datetime(2019, 1, 1)).days + 1)
		
		if val[5].isdigit() == False:
			try:
				n = "/aec/wf/2019/0"+day_of_year+"/"+str(val[5])+".*Z.20190"+day_of_year+"000000+"
				tr = obspy.read(n)
				
				tr[0].trim(tr[0].stats.starttime +(h *60 *60) + (mins * 60) + secs - tim, tr[0].stats.starttime +(h *60 *60) + (mins * 60) + secs + tim)
				data = tr[0][0:-1]
				fs = int(tr[0].stats.sampling_rate)					
				t                  = tr[0].times()
			except:
				continue
			
		else:	
			month2 = str(month1)
			if day1 < 10:
				day1 = '0'+str(day1)	
			if h < 23:			
				day2 = str(day1)
				if h < 10:
					h_u = '0'+str(h+1)
					h = '0'+str(h)
				else:
					h_u = str(h+1)
					h = str(h)
				
			else:
				h_u = '00'
				if month1 == '2' and day1 == '28':
					month2 = '03'
					day2 = '01'
				else:	
					if int(day1) >= 10:
						day2 = str(int(day1) + 1)
					else:
						day2 = '0'+str(int(day1)+1)
			try:		
				n = "/scratch/naalexeev/NODAL/2019-0"+str(month1)+"-"+str(day1)+"T"+str(h)+":00:00.000000Z.2019-0"+month2+"-"+day2+"T"+h_u+":00:00.000000Z."+station+".mseed"
				tr = obspy.read(n)
			
				tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
				data = tr[2][0:-1]
				fs = int(tr[2].stats.sampling_rate)				
				t                  = tr[2].times()
			except:
				continue
		base_dir = Path('input/plane_overtones/a'+plane[0:-1]+'.txt')
		
		try:
			# Time array
			t = np.arange(len(data)) / fs
			g = fs*240

			# Compute spectrogram
			frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9) 
			a, b = Sxx.shape
				
			MDF = np.zeros((a,b))
			for row in range(len(Sxx)):
				m = len(Sxx[row])
				p = sorted(Sxx[row])
				median = p[int(m/2)]
				for col in range(m):
					MDF[row][col] = median
			spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))
				
			# Find the index of the middle frequency
			middle_index = len(times) // 2

			# Extract the middle line of the spectrogram
			middle_column = spec[:, middle_index]
			peaks, _ = signal.find_peaks(middle_column, prominence=20) 
			np.diff(peaks)
			print(peaks, len(peaks))
			if os.path.exists(base_dir):
				output = open(base_dir,'a')
			else:
				output = open(base_dir,'w')
			for g in range(0,len(peaks)):
				output.write(str(station)+','+str(peaks[g])+','+str(10 * np.log10(middle_column[peaks[g]]))+'\n') #station,freq of peaks,relative amp
			output.close()
				
		except:
			continue

