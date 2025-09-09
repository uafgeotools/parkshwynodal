import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
import pytz
import obspy
import math
import matplotlib.dates as mdates
from pathlib import Path

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()

text = open('all_station_crossing_db.txt', 'r')

for line in text.readlines():
	val = line.split(',')	
	ht = datetime.datetime.utcfromtimestamp(int(val[2]))
	if ht.month == 2 and ht.day >= 11 or ht.month == 3:
		day_of_year = str((ht - datetime.datetime(2019, 1, 1)).days + 1)
	
		if val[5].isdigit() == False:
			try:
				n = "/aec/wf/2019/0"+day_of_year+"/"+str(val[5])+".*Z.20190"+day_of_year+"000000+"
				tr = obspy.read(n)

				h = ht.hour
				mins = ht.minute
				secs = ht.second

				tm = 300

				tr[0].trim(tr[0].stats.starttime + (h *60 *60) + (mins * 60) + secs - tm, tr[0].stats.starttime + (h *60 *60) +  (mins * 60) + secs + tm)

				t                  = tr[0].times()
				sampling_frequency = tr[0].stats.sampling_rate

				title    = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'

				fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     

				ax1.set(xlim=(0, 240), ylim=(-2000, 2000))
				ax1.plot(t, tr[0], 'k', linewidth=0.5)
				ax1.set_title(title)


				s,f,tm,im = ax2.specgram(tr[0], Fs=sampling_frequency, noverlap=int(0.8*256), cmap='hsv', detrend = 'linear', scale='dB', vmin = -120, vmax = 50)
				ax2.set_xlabel('Time - Seconds')

				ax2.set_ylabel('Frequency (Hz)')

				ax3 =fig.add_axes([0.88, 0.1, 0.02, 0.37])
				plt.colorbar(mappable=im, cax=ax3, spacing='uniform', label='Relative Amplitude (dB)')
				BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/perm_spec/'+ str(val[0]) + '/'+str(val[1])+'/'+str(val[5])+'/'
				make_base_dir(BASE_DIR)
				plt.savefig('/scratch/irseppi/nodal_data/plane_info/perm_spec/'+ str(val[0]) + '/'+str(val[1])+'/'+str(val[5])+'/'+str(val[1])+'_' + str(val[2]) + '.png')
				plt.close()
			except:
				print(day_of_year+' '+val[5])
				continue
	text.close()

