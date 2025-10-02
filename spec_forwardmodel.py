import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
import obspy
import datetime

flight_num = [530342801,528485724,528473220]
time = [1551066051,1550172833,1550168070]
sta = [1022,1272,1173]
day = [25,14,14]
month = [2,2,2]

forward_model = True

for n in range(0,3):
	if sta[n] == 'CCB' or sta[n] == 'F6TP' or sta[n] == 'F4TN' or sta[n] == 'F3TN' or sta[n] == 'F7TV':
		continue
	ht = datetime.datetime.utcfromtimestamp(time[n])
	mins = ht.minute
	secs = ht.second
	h = ht.hour
	tim = 120	
	h_u = h+1
	if h < 23:			
		day2 = day[n]
		if h < 10:
			h_u = '0'+str(h+1)
			h = '0'+str(h)
		else:
			h_u = h+1
			h = h
	else:
		h_u = '00'
		day2 = str(day[n]+1)
	if len(str(day[n])) == 1:
		day[n] = '0'+str(day[n])
		day2 = day[n]
	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/20190'+str(month[n])+str(day[n])+'_positions/20190'+str(month[n])+str(day[n])+'_'+str(flight_num[n])+'.csv', sep=",")

	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	tm = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	head = flight_data['heading']
	for line in range(len(tm)):
		if str(tm[line]) == str(time[n]):
			speed_knots = flight_data['speed'][line]
			speed_mps = speed_knots * 0.514444
			alt_ft = flight_data['altitude'][line]
			alt_m = alt_ft * 0.3048

	p = "/scratch/naalexeev/NODAL/2019-0"+str(month[n])+"-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-0"+str(month[n])+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta[n])+".mseed"
	tr = obspy.read(p)
	tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
	data = tr[2][0:-1]
	fs = int(tr[2].stats.sampling_rate)
	title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
	t = tr[2].times()
	# Time array
	t = np.arange(len(data)) / fs
	g = fs*240
	# Compute spectrogram
	frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 

	a, b = Sxx.shape

	MDF = np.zeros((a,b))
	for row in range(len(Sxx)):
		m = len(Sxx[row])
		p = sorted(Sxx[row])
		median = p[int(m/2)]
		for col in range(m):
			MDF[row][col] = median 
	
	start_time = (tr[2].stats.starttime + (mins * 60) + secs - tim).timestamp
	print('Start time:', start_time)
	
	spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

	# Find the index of the middle frequency
	middle_index = len(times) // 2
	middle_column = spec[:, middle_index]
	vmin = 0  
	vmax = np.max(middle_column) 

	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

	ax1.plot(t, data, 'k', linewidth=0.5)
	ax1.set_title(title)
	
	ax1.margins(x=0)

	if n == 0:
		fnot = [115]
		tpr = np.arange(0, 241, 1)
		c = 320
		v0 = 61.22 
		tprime0 = (((1551066047.29- (time[n] - 120)))) + np.sqrt(1655.578**2 + 426.72**2)/c -2 
		l = np.sqrt(1655.58**2 + 426.72**2) 

	if n == 1:
		fnot = [115]
		tpr = np.arange(0, 241, 1)
		c = 320
		v0 = 95.429362
		tprime0 = np.abs(((time[n] - 120) -  1550172810.90)) + np.sqrt(521.55**2 + 2933.7**2)/c 
		l = np.sqrt(521.55**2 + 2933.7**2)

	if n == 2:
		fnot = [131]
		tprime0 = 100
		tpr = np.arange(0, 241, 1)
		c = 320
		v0 = 123.46656
		l = np.sqrt(660.03**2 + 4876.8**2)

	print(l)
	print(c)

	tprime0 = l/c
	t = ((tprime0 - tprime0)+(l/c)- np.sqrt((v0/c)**2*(tprime0-tprime0)**2 + (v0/c)**2*(2*l/c)*((tprime0-tprime0)+l**2/c**2)))/(1-v0**2/c**2)
	print('tprime0:', tprime0)
	print('t:', t)


	for f0 in fnot:
		ft = []
		for tprime in tpr:
			t = ((tprime - tprime0)+(l/c)- np.sqrt((v0/c)**2*(tprime-tprime0)**2 + (v0/c)**2*(2*l/c)*(tprime-tprime0)+l**2/c**2))/(1-v0**2/c**2)
			ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))								
			ft.append(ft0p)

		ax2.plot(tpr, ft, 'g', linewidth=0.5)
		ax2.set_title("Forward Model: t'= "+str(tprime0)+' sec, v0 = '+str(v0)+' m/s, l = '+str(l)+' m, \n' + 'f0 = '+str(fnot)+' Hz', fontsize='x-small')
	# Plot spectrogram
	cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
	ax2.set_xlabel('Time (s)')


	ax2.axvline(x=tprime0, c = 'g', ls = '--', label='Estimated arrival: '+str(tprime0)+' s')
	ax2.legend(loc='upper right',fontsize = 'x-small')
	ax2.set_ylabel('Frequency (Hz)')

	
	ax2.margins(x=0)
	ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

	plt.colorbar(mappable=cax, cax=ax3)
	ax3.set_ylabel('Relative Amplitude (dB)')

	ax2.margins(x=0)
	ax2.set_xlim(0, 240)
	
	# Plot overlay
	spec2 = 10 * np.log10(MDF)
	middle_column2 = spec2[:, middle_index]
	vmin = np.min(middle_column2)
	vmax = np.max(middle_column2)

	# Create ax4 and plot on the same y-axis as ax2
	ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) 
	ax4.plot(middle_column2, frequencies, c='orange')  
	ax4.set_ylim(0, int(fs/2))
	ax4.set_xlim(vmax*1.1, vmin) 
	ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
	ax4.grid(axis='y')
	
	plt.show()


