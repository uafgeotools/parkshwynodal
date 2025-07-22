import numpy as np
import numpy.linalg as la
import obspy
from matplotlib import pyplot as plt
from datetime import datetime, timezone
from prelude import calc_ft, S, speed_of_sound
from scipy.signal import spectrogram
from plot_func import  remove_median

####################################################################################################################################################################################################################################################################################################################

def df(f0,v0,l,tp0,tp,c):   
    """
	Calculate the derivatives of f with respect to f0, v0, l, and tp0.

	Parameters:
	f0 (float): Fundamental frequency produced by the aircraft.
	v0 (float): Velocity of the aircraft.
	l (float): Distance of closest approach between the station and the aircraft.
	tp0 (float): Time of that the central frequency of the overtones occur, when the aircraft is at the closest approach to the station.
	tp (float): Array of times.
	Returns:
	tuple: A tuple containing the derivatives of f with respect to f0, v0, l, and tp0.
	"""
    #print('f0: ', f0, 'v0: ', v0, 'l: ', l, 'tp0: ', tp0)
    #derivative with respect to f0
    f_derivef0 = (1 / (1 - (c * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2  + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))) /((c**2 - v0**2) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2))))

    #derivative of f with respect to v0
    f_derivev0 = -(f0 * v0 * (-2 * l**4 * v0**4 + l**2 * (tp - tp0)**2 * v0**6 + c**6 * (tp - tp0) * (2 * l**2 + (tp - tp0)**2 * v0**2) * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) + 
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
    f_derivetprime0 = ((f0 * l**2 * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) / 
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

    return f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec

#####################################################################################################################################################################################################################################################################################################################

def invert_f(mprior, coords_array, num_iterations,sigma = 10):
	"""
	Inverts the function f using the given initial parameters and data array.

	Args:
		m0 (numpy.ndarray): Initial parameters for the function f.
		coords_array (numpy.ndarray): Data picks along overtone doppler curve.
		num_iterations (int): Number of iterations to perform.

	Returns:
		numpy.ndarray: The inverted parameters for the function f.
	"""
	off_diagonal = False #if True, the covariance matrix will have off-diagonal elements, if False, it will be diagonal
	dw,_ = coords_array.shape
	fobs = coords_array[:,1]
	tobs = coords_array[:,0]
	n = 0

	f0_prior = 80
	v0_prior = 100
	l_prior = 10000
	tprime0_prior = 40
	c_prior = 60
	cprior0 = np.zeros((5,5))

	cprior0[0][0] = f0_prior**2
	cprior0[1][1] = v0_prior**2
	cprior0[2][2] = l_prior**2
	cprior0[3][3] = tprime0_prior**2
	cprior0[4][4] = c_prior**2
	if off_diagonal:
		cprior0[0][3] =  -0.4*f0_prior*tprime0_prior

		cprior0[1][2] = -0.7*v0_prior*l_prior
		cprior0[1][4] = 0.85*v0_prior*c_prior
		
		cprior0[2][1] = -0.7*v0_prior*l_prior
		cprior0[2][4] = -0.7*l_prior*c_prior

		cprior0[3][0] =  -0.4*f0_prior*tprime0_prior
		
		cprior0[4][1] = 0.85*v0_prior*c_prior
		cprior0[4][2] = -0.7*l_prior*c_prior     

	cprior = cprior0 * (5)

	Cd0 = np.zeros((len(fobs), len(fobs)), int)
	np.fill_diagonal(Cd0, sigma**2)
	Cd = Cd0*(dw)
	mnew = mprior.copy() #mprior is the initial guess for the parameters, mnew is the updated guess
	while n < num_iterations:
		m = mnew
		fpred = []
		G = np.zeros((dw,5)) #partial derivative matrix of f with respect to m
		#partial derivative matrix of f with respect to m 
		for i in range(0,dw):
			f0 = m[0]
			v0 = m[1]
			l = m[2]
			tprime0 = m[3]
			c = m[4]
			tprime = tobs[i]
			t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
			ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
			f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec = df(m[0], m[1], m[2], m[3], tobs[i],m[4])
			
			G[i,0:5] = [f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec]

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

		n += 1
		print(mnew)
	covmlsq = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior0))
	F_m = S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
	return mnew, covmlsq, F_m

#################################################################################################################################################################################################################################################################################################################################################################

c = speed_of_sound(-33)

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

#insert method to get initial model here
f0 = coords_array[1,1]+coords_array[2,1]/2
tprime0 = coords_array[0,0] 
v0 = c*abs(coords_array[1,1]-coords_array[2,1]) / (2 * f0)
slope = (coords_array[4,1] - coords_array[3,1]) / (coords_array[4,0] - coords_array[3,0])
l = -((f0*v0**2/c)*(1-(v0/c)**2)**(-3/2))/slope #(c**2*f0*v0**2*np.sqrt(c**2 - v0**2)/(c**2 - v0**2)**2)/abs(slope)

m0 = [f0, v0, l, tprime0,c]
print('Initial model:', m0)
m, covm, F_m = invert_f(m0, coords_array, num_iterations=8)
ft = calc_ft(times, m[3], m[0], m[1], m[2], m[4])

peaks = []
coord_inv = []
corridor_width = 6 
for t_f in range(len(times)):
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

coord_inv_array = np.array(coord_inv)

m,_,F_m = invert_f(m, coord_inv_array, num_iterations=4)

ft = calc_ft(times, m[3], m[0], m[1], m[2], m[4])

delf = np.array(ft) - np.array(peaks)

new_coord_inv_array = []
for i in range(len(delf)):
    if np.abs(delf[i]) <= 3:
        new_coord_inv_array.append(coord_inv_array[i])
coord_inv_array = np.array(new_coord_inv_array)

m,covm,F_m = invert_f(m0, coord_inv_array, num_iterations=20, sigma=3)

f0 = m[0]
v0 = m[1]
l = m[2]
tprime0 = m[3]
c = m[4]

nx,ny = covm.shape


sigma = np.sqrt(np.diag(covm))
outer_v = np.outer(sigma,sigma)
Crho = covm / outer_v

Crho[covm == 0] = 0
'''
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
'''
covm = np.sqrt(np.diag(covm))
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

