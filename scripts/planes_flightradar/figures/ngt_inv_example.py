import numpy as np
import numpy.linalg as la
import obspy
from matplotlib import pyplot as plt
from datetime import datetime, timezone
from src.doppler_funcs import calc_ft, invert_f, speed_of_sound
from scipy.signal import spectrogram
from src.main_inv_fig_functions import  remove_median

generate_samples = False
c = speed_of_sound(-33.7)
start_time = 1550158642.26246    
ht = datetime.fromtimestamp(start_time, tz=timezone.utc)                      
h = ht.hour
mins = ht.minute
secs = ht.second
month = ht.month
day = ht.day
h_u = str(h+1)

p = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day)+"T"+str(h_u)+":00:00.000000Z.1173.mseed"
tr = obspy.read(p)
tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs , tr[2].stats.starttime + (mins * 60) + secs + 240)
data = tr[2][:]
fs = int(tr[2].stats.sampling_rate)
title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
torg = tr[2].times()

# Compute spectrogram
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 

spec, MDF = remove_median(Sxx)

middle_index =  len(times) // 2
middle_column = spec[:, middle_index]
vmin = 0  
vmax = np.max(middle_column) 

x = [112.48911983478979, 59.65932080234049, 186.52395930932946, 102.98341205040444, 120.34960896418536]
y = [140.02964002964, 188.29218829218826, 93.7170937170937, 153.9234039234039, 128.81712881712878]

coords = [(x[i], y[i]) for i in range(len(x))]
coords_array = np.array(coords)
if generate_samples:
	fig_num = 6
else:
	fig_num = 5
# Create a subplot for the visualization
fig, ax = plt.subplots(fig_num,1,figsize=(8/1.4, 14/1.4),sharex=False)
cax = ax[0].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[0].axhline(y=188.29218829218826, color='black', linestyle='--', linewidth=1)
ax[0].axhline(y=93.7170937170937, color='black', linestyle='--', linewidth=1)
ax[0].axhline(y=(188.29218829218826+93.7170937170937)/2, color='red', linestyle='--', linewidth=0.7)
ax[0].axvline(x=112.48911983478979, color='red', linestyle='--', linewidth=0.7)
slope = (coords_array[4,1] - coords_array[3,1]) / (coords_array[4,0] - coords_array[3,0])
#Create dashed line with the slope that is tangent to the first point
x_values = [coords_array[3,0], coords_array[4,0]]
y_values = [coords_array[3,1], coords_array[4,1]]
# Plot the dashed line before scatter plots so it appears underneath
ax[0].plot(x_values, y_values, color='blue', linestyle='--', linewidth=1, zorder=1)

# Add points at x=70 and x=150 using the slope
y_70 = coords_array[3,1] + slope * (70 - coords_array[3,0])
y_150 = coords_array[3,1] + slope * (150 - coords_array[3,0])
ax[0].plot([70, 150], [y_70, y_150], color='blue', linestyle='--', linewidth=1, zorder=1)
# Move scatter plots after all lines so they appear on top
ax[0].scatter(coords_array[1:3, 0], coords_array[1:3, 1], c='black', marker='x', s=100, linewidths=3,label="f_initial + f_final")
ax[0].scatter(coords_array[0, 0], coords_array[0, 1], c='red', marker='x', s=100, linewidths=3, label="t'0 + f0")
ax[0].scatter(coords_array[3:5, 0], coords_array[3:5, 1], c='blue', marker='x', s=100, linewidths=3, label="Slope of l")
ax[0].set_ylabel('Frequency (Hz)')
ax[0].set_title("(a) data picks with steps to get prior model", fontsize='small')

cax = ax[1].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

#insert method to get initial model here
f0 = (coords_array[1,1]+coords_array[2,1])/2 - 20
tprime0 = coords_array[0,0] 
v0 = c*abs(coords_array[1,1]-coords_array[2,1]) / (2 * f0)
slope = (coords_array[4,1] - coords_array[3,1]) / (coords_array[4,0] - coords_array[3,0])
l = -((f0*v0**2/c)*(1-(v0/c)**2)**(-3/2))/slope #(c**2*f0*v0**2*np.sqrt(c**2 - v0**2)/(c**2 - v0**2)**2)/abs(slope)

m0 = [f0, v0, l, tprime0,c]
prior_sigma = [50, 70, 2000, 30, 60] #initial prior sigma values for f0, v0, l, tprime0, c


print('Initial model:', m0)
ft = calc_ft(times, m0[3], m0[0], m0[1], m0[2], m0[4])
ax[1].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[1].scatter(coords_array[1:3, 0], coords_array[1:3, 1], c='black', marker='x', s=100, linewidths=3,label="f_initial + f_final")
ax[1].scatter(coords_array[0, 0], coords_array[0, 1], c='red', marker='x', s=100, linewidths=3, label="t'0 + f0")
ax[1].scatter(coords_array[3:5, 0], coords_array[3:5, 1], c='blue', marker='x', s=100, linewidths=3, label="Slope of l")
ax[1].set_ylabel('Frequency (Hz)')
ax[1].set_title("(b) initial prior model", fontsize='small')
m, covm,_, F_m = invert_f(m0, prior_sigma, coords_array, num_iterations=5)
ft = calc_ft(times, m[3], m[0], m[1], m[2], m[4])

cax = ax[2].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[2].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[2].set_ylabel('Frequency (Hz)')
ax[2].set_title("(c) inverted model => updated prior model", fontsize='small')
peaks = []
coord_inv = []
upper_array = []
lower_array = []
corridor_width = 10
time_corr = np.arange(0, 240, 1)
for ttt in time_corr:
    t_f = (np.abs(times - ttt)).argmin()
    upper = int(ft[t_f] + corridor_width)
    lower = int(ft[t_f] - corridor_width)
    if lower < 0:
        lower = 0
    if upper > len(frequencies):
        upper = len(frequencies)
    tt = spec[lower:upper, t_f]
    try:
        max_amplitude_index = np.argmax(tt)
    except:
        continue
    max_amplitude_frequency = frequencies[max_amplitude_index+lower]
    peaks.append(max_amplitude_frequency)
    coord_inv.append((times[t_f], max_amplitude_frequency))
    upper_array.append(upper)
    lower_array.append(lower)

coord_inv_array = np.array(coord_inv)
cax = ax[3].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[3].plot(coord_inv_array[:, 0], np.array(upper_array), 'r', linewidth=1)
ax[3].plot(coord_inv_array[:, 0], np.array(lower_array), 'r', linewidth=1)
ax[3].set_title("(d) data extracted from model corridor (updated prior model \u00B1 10)", fontsize='small')

prior_sigma = [5,10,600,5,30] #prior sigma values for f0, v0, l, tprime0, c
m,_,_,F_m = invert_f(m, prior_sigma, coord_inv_array, num_iterations=3)

f0 = m[0]
v0 = m[1]
l = m[2]
tprime0 = m[3]
c = m[4]
print(l)
l = -((f0*v0**2/c)*(1-(v0/c)**2)**(-3/2))/slope
print(l)
m[2] = l
ft = calc_ft(time_corr, m[3], m[0], m[1], m[2], m[4])

delf = np.array(ft) - np.array(peaks)

new_coord_inv_array = []
for i in range(len(delf)):
    if np.abs(delf[i]) <= 3:
        new_coord_inv_array.append(coord_inv_array[i])
coord_inv_array = np.array(new_coord_inv_array)

ax[3].scatter(coord_inv_array[:, 0], coord_inv_array[:, 1], c='black', marker='x', s=20)
ax[3].set_ylabel('Frequency (Hz)')

prior_sigma = [5,10,500,30,80] #prior sigma values for f0, v0, l, tprime0, c
m,covm0,covm_norm,F_m = invert_f(m, prior_sigma, coord_inv_array, num_iterations=6, sigma=3)

f0 = m[0]
v0 = m[1]
l = m[2]
tprime0 = m[3]
c = m[4]

ft = calc_ft(times, tprime0, f0, v0, l, c)
cax = ax[4].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax[4].plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=1) 
ax[4].set_ylabel('Frequency (Hz)')
ax[4].set_title("(e) posterior model", fontsize='small')

if generate_samples:
	nx,ny = covm0.shape
	# initialize samples
	covm_samples = np.empty((5,1000))
	m_samples = np.zeros((5,1000))
	ft_matrix = np.zeros((1000,len(times)))
	# generate samples of the posterior
	R = np.linalg.cholesky(covm0)
	cax = ax[5].pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
	for jj in range(1000):
		covm_samples[:,jj] = (R @ np.random.randn(5,1)).flatten()
		m_samples[:,jj] = covm_samples[:,jj] + m.flatten()
		f0_s = m_samples[0,jj]
		v0_s = m_samples[1,jj]
		l_s = m_samples[2,jj]
		tprime0_s = m_samples[3,jj]
		c_s = m_samples[4,jj]
		ft = calc_ft(times, tprime0_s, f0_s, v0_s, l_s, c_s)
		ft_matrix[jj, :] = ft 
	ax[5].set_ylim(0, 250)
	std_samples  = np.std(ft_matrix,axis=0)
	ft = calc_ft(times, tprime0, f0, v0, l, c)
	ax[5].plot(times, ft+std_samples, color='red', linewidth=0.5)
	ax[5].plot(times, ft-std_samples, color='red', linewidth=0.5)
	ax[5].set_title("(f) standard deviation of posterior model samples", fontsize='small')

#make all axis tick labels smaller
for i in range(fig_num):
	ax[i].tick_params(axis='both', which='major', labelsize='x-small')
	ax[i].tick_params(axis='both', which='minor', labelsize='x-small')
#make  the gap between subplots smaller
plt.subplots_adjust(hspace=0.3)
ax[fig_num-1].set_xlabel('Time (s)')
plt.tight_layout()
fig.savefig("inversion_steps.jpg", dpi=600)
plt.close()



make_final_plot = False
if make_final_plot:
	covm = np.sqrt(np.diag(covm0))
	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

	ax1.plot(torg, data, 'k', linewidth=0.5)
	ax1.set_title(title)

	ax1.margins(x=0)
	ax1.set_position([0.125, 0.6, 0.775, 0.3])  # Move ax1 plot upwards

	# Plot spectrogram
	cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
	ax2.set_xlabel('Time (s)')

	ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label= "t\u2080' = " + "%.2f" % tprime0 +' s')

	ft = calc_ft(times, tprime0, f0, v0, l, c)

	ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 
	tprime = tprime0
	t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
	ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

	ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 

	fss = 'x-small'

	ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % covm[3] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % covm[1] + ' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % covm[4] + ' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % covm[2] + ' m, \n' + 'f\u2080 =' + "%.2f" % f0 + ' \u00B1 ' + "%.2f" % covm[0] +' Hz\nMisfit: ' + "%.4f" % F_m, fontsize=fss)


	ax2.legend(loc='upper right',fontsize = 'small')
	ax2.set_ylabel('Frequency (Hz)')

	ax2.margins(x=0)
	ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

	plt.colorbar(mappable=cax, cax=ax3)
	ax3.set_ylabel('Relative Amplitude (dB)')

	ax2.margins(x=0)
	ax2.set_xlim(0, 240)
	ax2.set_ylim(0, int(fs/2))

	# Plot overlay
	spec2 = 10 * np.log10(MDF)
	middle_column2 = spec2[:, middle_index]
	vmin2 = np.min(middle_column2)
	vmax2 = np.max(middle_column2)

	# Create ax4 and plot on the same y-axis as ax2
	ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) 
	ax4.plot(middle_column2, frequencies, c='#ff7f00')  
	ax4.set_ylim(0, int(fs/2))
	ax4.set_xlim(vmax2*1.1, vmin2) 
	ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
	ax4.grid(axis='y')
	plt.show()
	plt.close()

	if generate_samples:
		plt.figure()
		plt.plot(times, 2*std_samples)
		plt.ylim(0,2)
		plt.show()
		plt.close()

	plot_posterior = False
	if plot_posterior:
		sigma = np.sqrt(np.diag(covm_norm))
		outer_v = np.outer(sigma,sigma)
		Crho = covm_norm / outer_v

		Crho[covm_norm == 0] = 0

		gridlines=False
		colormap='seismic'
		plt.figure(figsize=(10, 10))
		plt.imshow(Crho,cmap=colormap)
		plt.xticks(ticks=range(np.shape(Crho)[1]),labels=[str(val) for val in range(1,np.shape(Crho)[1]+1)])
		plt.yticks(ticks=range(np.shape(Crho)[0]),labels=[str(val) for val in range(1,np.shape(Crho)[0]+1)])
		if gridlines:
			xgrid = np.array(range(np.shape(Crho)[1] + 1)) - 0.5
			ygrid = np.array(range(np.shape(Crho)[0] + 1)) - 0.5
			for gridline in xgrid:
				plt.axvline(x=gridline,color='k',linewidth=1)
			for gridline in ygrid:
				plt.axhline(y=gridline,color='k',linewidth=1)
		plt.colorbar()
		plt.show()
		plt.close()
