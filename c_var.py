import numpy as np
import pandas as pd
import json
import obspy
import datetime
from datetime import datetime, timezone
from pyproj import Proj
from prelude import *
from scipy.signal import find_peaks, spectrogram
from plot_func import *
from obspy.clients.nrl import NRL
import os
##############################################################################################################################################################################################################

def plot_spectrgram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, covm, flight, middle_index, tarrive, closest_time, dir_name, plot_show=True):
    """
    Plot the spectrogram and spectrum of the given data.

    Args:
        data (array): The waveform data.
        fs (int): The sampling frequency.
        torg (array): The time array.
        title (str): The title of the plot.
        spec (array): The spectrogram data.
        times (array): The time array for the spectrogram.
        frequencies (array): The frequency array for the spectrogram.
        tprime0 (float): The estimated arrival time.
        v0 (float): The velocity.
        l (float): The distance.
        c (float): The speed of sound.
        f0_array (array): The array of frequencies.
        arrive_time (array): The arrival time array.
        MDF (array): Median removed from spectrogram.
        covm (array): The covariance matrix.
        flight (int): The flight number.
        middle_index (int): The index of the middle column.
        closest_time (float): The closest time.
        dir_name (str): The directory name.

    Returns:
        str: The user assigned quality number.
    """
    # Plot settings and calculations
    vmin = np.min(arrive_time) #tprime0
    vmax = np.max(arrive_time)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

    ax1.plot(torg, data, 'k', linewidth=0.5)
    ax1.set_title(title)

    ax1.margins(x=0)
    ax1.set_position([0.125, 0.6, 0.775, 0.3])  # Move ax1 plot upwards

    # Plot spectrogram
    cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
    ax2.set_xlabel('Time (s)')
    f0lab = []
    #ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label= "t\u2080' = " + str(np.round(tprime0,2))+' s')
    ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label= "t\u2080' = " + "%.2f" % tprime0 +' s')
    for pp in range(len(f0_array)):
        f0 = f0_array[pp]
        
        ft = calc_ft(times, tprime0, f0, v0, l, c)

        ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 
        tprime = tprime0
        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
        
        ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 

    fss = 'x-small'
    f0lab = sorted(f0_array)

    if len(f0_array) <= 1:
        med_df = "NaN"
        mad_df = "NaN"

    else:
        #Generate random samples of f0 values withing their sigma from the covariance matrix 
        #Calculate the median of the differences and MAD to obtain error
        f_range = []
        NTRY = 1000
        for N in range(NTRY):
            ftry = []
            for c_index  in range(3, len(covm)):
                xmin = f0_array[c_index-4] - covm[c_index]
                xmax = f0_array[c_index-4] + covm[c_index]
                xtry = xmin + (xmax-xmin)*np.random.rand()
                ftry.append(xtry)

            ftry = np.sort(ftry)
            f1 = []
            for g in range(len(ftry)):
                if g == 0:
                    continue
                diff = ftry[g] - ftry[g - 1]
                f1.append(diff)
            med = np.nanmedian(f1)
            f_range.append(med)
        med_df = np.nanmedian(f_range)
        mad_df = np.nanmedian(np.abs(f_range - med_df))
    if med_df == "NaN":
        ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % covm[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % covm[0]+' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % covm[3]+' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % covm[1] + ' m, \n' + 'f\u2080 = ['+', '.join(["%.2f" % f for f in f0lab]) +'] \u00B1 ' + "%.2f" % np.median(covm[3:]) +' Hz, df\u2080 = NaN \u00B1 NaN Hz\nMisfit: ' + "%.4f" % F_m, fontsize=fss)
    else:
        ax2.set_title("t\u2080'= "+ "%.2f" % tprime0 + ' \u00B1 ' + "%.2f" % covm[2] + ' s, v\u2080 = ' + "%.2f" % v0 +' \u00B1 ' + "%.2f" % covm[0] + ' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % covm[3] + ' m/s, l = '+ "%.2f" % l +' \u00B1 ' + "%.2f" % covm[1] + ' m, \n' + 'f\u2080 = ['+', '.join(["%.2f" % f for f in f0lab]) +'] \u00B1 ' + "%.2f" % np.median(covm[3:]) +' Hz, df\u2080 = ' + "%.2f" % med_df + ' \u00B1 ' + "%.2f" % mad_df + ' Hz\nMisfit: ' + "%.4f" % F_m + ' [FH/VT]', fontsize=fss)
    ax2.axvline(x=tarrive, c = '#e41a1c', ls = '--',linewidth=0.5,label= r'$t_{i}$ = ' + "%.2f" % tarrive +' s')

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
    if plot_show:
        plt.show()     
        qnum = input('What quality number would you give this?(first num for data quality(0-3), second for ability to fit model to data(0-1))')
    else:
        qnum = '__'
    
    fig.savefig(dir_name+'/'+str(closest_time)+'_'+str(flight)+'.png')
    plt.close()

    return qnum

##############################################################################################################################################################################################################

def plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, dir_name):
    """
    Plot the spectrum with arrival time markers.

    Args:
        spec (numpy.ndarray): The spectrum data.
        frequencies (numpy.ndarray): The frequencies.
        tprime0 (float): The reference arrival time.
        v0 (float): The velocity.
        l (float): The distance.
        c (float): The speed of light.
        f0_array (list): The list of frequencies.
        arrive_time (numpy.ndarray): The arrival times.
        fs (int): The sampling frequency.
        closest_index (int): The index of the closest time.
        closest_time (float): The closest time.
        sta (int or str): The station identifier.
        dir_name (str): The directory name.

    Returns:
        None
    """
    vmax = np.max(arrive_time)
    fig = plt.figure(figsize=(10,6))
    plt.grid()

    plt.plot(frequencies, spec[:,closest_index], c='#377eb8')
        
    for pp in range(len(f0_array)):
        f0 = f0_array[pp]
        if fs/2 < f0:
            continue
        tprime = tprime0
        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
        if ft0p == np.nan:
            continue
        if ft0p > 250:
            continue
        try:
            upper = int(ft0p + 10)
            lower = int(ft0p - 10)
            tt = spec[lower:upper, closest_index]
            if upper > 250:
                freqp = ft0p
                ampp = np.interp(ft0p, frequencies, arrive_time)
            elif lower < 0:
                freqp = ft0p
                ampp = np.interp(ft0p, frequencies, arrive_time)
            else:
                ampp = np.max(tt)
                freqp = np.argmax(tt)+lower
            plt.scatter(freqp, ampp, color='black', marker='x', s=100)
            if isinstance(sta, int):
                plt.text(freqp - 10, ampp + 0.8, freqp, fontsize=17, fontweight='bold')
            else:
                plt.text(freqp - 1, ampp + 0.8, freqp, fontsize=17, fontweight='bold')  
        except:
            continue
    plt.xlim(0, int(fs/2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0,vmax*1.1)
    plt.xlabel('Frequency (Hz)', fontsize=17)
    plt.ylabel('Relative Amplitude at t = {:.2f} s (dB)'.format(tprime0), fontsize=17)

    fig.savefig(dir_name + '/'+str(sta)+'_' + str(closest_time) + '.png')
    plt.close()
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

def invert_f(m0, coords_array, num_iterations,sigma = 10):
	"""
	Inverts the function f using the given initial parameters and data array.

	Args:
		m0 (numpy.ndarray): Initial parameters for the function f.
		coords_array (numpy.ndarray): Data picks along overtone doppler curve.
		num_iterations (int): Number of iterations to perform.

	Returns:
		numpy.ndarray: The inverted parameters for the function f.
	"""
	w,_ = coords_array.shape
	fobs = coords_array[:,1]
	tobs = coords_array[:,0]
	m = m0
	n = 0
	while n < num_iterations:
		fnew = []
		G = np.zeros((w,5)) #partial derivative matrix of f with respect to m
		#partial derivative matrix of f with respect to m 
		for i in range(0,w):
			f0 = m[0]
			v0 = m[1]
			l = m[2]
			tprime0 = m[3]
			c = m[4]
			tprime = tobs[i]
			t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
			ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
			f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec = df(m[0], m[1], m[2], m[3], tobs[i],c)
			
			G[i,0:5] = [f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec]

			fnew.append(ft0p) 
		try:
			covmlsq = (sigma**2)*la.inv(G.T@G)
		except:
			covmlsq = (sigma**2)*la.pinv(G.T@G)
		try:
			m = np.reshape(np.reshape(m0,(5,1))+ np.reshape(la.inv(G.T@G)@G.T@(np.reshape(fobs, (len(coords_array), 1)) - np.reshape(np.array(fnew), (len(coords_array), 1))), (5,1)), (5,))
		except:
			m = np.reshape(np.reshape(m0,(5,1))+ np.reshape(la.pinv(G.T@G)@G.T@(np.reshape(fobs, (len(coords_array), 1)) - np.reshape(np.array(fnew), (len(coords_array), 1))), (5,1)), (5,))
		print(m)
		m0 = m
		n += 1
	F_m = Sd(fnew, fobs, len(fobs), sigma)
	return m, covmlsq, F_m

#####################################################################################################################################################################################################################################################################################################################

def full_inversion(fobs, tobs, freqpeak, peaks, peaks_assos, tprime, tprime0, ft0p, v0, l, f0_array, mprior, c, w, num_iterations = 4, sigma = 10):
	"""
	Performs inversion using all picked overtones. 

	Args:
		fobs (numpy.ndarray): Picked frequency values from individual overtone inversion picks.
		tobs (numpy.ndarray): Picked time values from individual overtone inversion picks.
		freqpeak (numpy.ndarray): Center time value of doppler curve for each overtone
		peaks (numpy.ndarray): Value of the frequency at the center time of the doppler curve for each overtone.

	Returns:
		numpy.ndarray: The inverted parameters for the function f. Velocity of the aircraft, distance of closest approach, time of closest approach, and the fundamental frequency produced by the aircraft.
		numpy.ndarray: The covariance matrix of the inverted parameters.
		numpy.ndarray: The array of the fundamental frequency produced by the aircraft.
	"""

	qv = 0

	cprior = np.zeros((w+4,w+4))

	for row in range(len(cprior)):
		if row == 0:
			cprior[row][row] = 20**2
		elif row == 1:
			cprior[row][row] = 800**2
		elif row == 2:
			cprior[row][row] = 80**2
		elif row == 3:
			cprior[row][row] = 10**2
		else:
			cprior[row][row] = 10**2
	
	Cd = np.zeros((len(fobs), len(fobs)), int)
	np.fill_diagonal(Cd, sigma**2)
	mnew = np.array(mprior)
	
	while qv < num_iterations:
		G = np.zeros((0,w+4))
		fnew = []
		cum = 0
		for p in range(w):
			new_row = np.zeros(w+4)
			f0 = f0_array[p]
			
			for j in range(cum,cum+peaks_assos[p]):
				tprime = tobs[j]
				t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
				ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

				f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec = df(f0,v0,l,tprime0, tobs[j],c)
			
				new_row[0] = f_derivev0
				new_row[1] = f_derivel
				new_row[2] = f_derivetprime0
				new_row[3] = f_derivec
				new_row[4+p] = f_derivef0
						
				G = np.vstack((G, new_row))
						
				fnew.append(ft0p)
		
			cum = cum + peaks_assos[p]

		m = np.array(mnew) + cprior@G.T@la.inv(G@cprior@G.T+Cd)@(np.array(fobs)- np.array(fnew))
		mnew = m
		v0 = mnew[0]
		l = mnew[1]
		tprime0 = mnew[2]
		c = mnew[4]
		f0_array = mnew[4:]

		print(m)
		qv += 1
	covm = la.inv(G.T@la.inv(Cd)@G + la.inv(cprior))
	F_m = Sd(fnew, fobs, len(fobs), sigma)
	return m, covm, f0_array, F_m


nrl = NRL()

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

rerun_fig = True #False #Flag rerun the figures without saving the inversion results = True

# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

for li in file_in.readlines():
    text = li.split(',')
    date = text[0]
    flight_num = text[1]
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate
    dist_m = float(text[4])   # Distance in meters
    closest_time = float(text[5])
    alt = float(text[6]) 
    speed_mps = float(text[7])  # Speed in meters per second
    sta = text[9]
    equip = text[10]

    if equip == 'B737' or equip == 'DH8A' or equip == 'AT73' or equip == 'B763':
        continue
    
    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    spec_dir = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
    
    if os.path.exists(spec_dir): 
        go = True
    else:
        continue

    for i in range(len(stations)):
        if stations[i] == sta:
            seismo_lat = seismo_latitudes[i]
            seismo_lon = seismo_longitudes[i]
            elev = elevations[i]
            break
    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)

    #if rerun_fig == False:
    output = open('output/with_c/' + equip + 'data_atmosphere_full.csv', 'a')

    input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(closest_time) + '_' + str(lat) + '_' + str(lon) + '.dat'
    
    try:
        file =  open(input_files, 'r') 
    except:
        print('No tempurature file for: ', date, flight_num, sta)
        continue

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
            for i in range(len(item['values'])):
                if abs(float(item['values'][i]) - float(alt/1000)) < hold:
                    hold = abs(float(item['values'][i]) - float(alt/1000))
                    z_index = i

    for item in data_list:
        if item['parameter'] == 'T':
            Tc = - 273.15 + float(item['values'][z_index])

    c = speed_of_sound(Tc)
    sound_speed = c
    tarrive = calc_time(closest_time,dist_m,alt,c) 

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    time = flight_data['snapshot_id']
    timestamps = flight_data['snapshot_id']
    speed = flight_data['speed']
    altitude = flight_data['altitude']

    #Must use the tarrive time to get the correct data
    ht = datetime.fromtimestamp(tarrive, tz=timezone.utc)
    if equip[0] == 'B' and equip[0:1] != 'BE':
        wind = 120
    else:
        wind = 120
    #start_time = tarrive - wind
    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'   
    if Path(file_name).exists():
        coords = []
        if Path(file_name).is_dir():
            continue
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')

                try:
                    start_time = float(pick_data[2])
                    print('Start time: ', start_time)
                    pp = 'yep'
                    break
                except:
                    pp = 'nope'
                    break
    ht = datetime.fromtimestamp(start_time+120, tz=timezone.utc)                        
    if pp == 'nope':
        continue
    h = ht.hour
    mins = ht.minute
    secs = ht.second
    month = ht.month
    day = ht.day

    h_u = str(h+1)
    if h < 23:			
        day2 = str(day)
        if h < 10:
            h_u = '0'+str(h+1)
            h = '0'+str(h)
        else:
            h_u = str(h+1)
            h = str(h)
    else:
        h_u = '00'
        day2 = str(day+1)
    if len(str(day)) == 1:
        day = '0'+str(day)
        day2 = day

    try:
        p = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta)+".mseed"
        tr = obspy.read(p)
        tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - wind, tr[2].stats.starttime + (mins * 60) + secs + wind)
        data = tr[2][:]
        fs = int(tr[2].stats.sampling_rate)
        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
        torg = tr[2].times()
        if len(data) == 0:
            data = tr[1][:]
            fs = int(tr[1].stats.sampling_rate)
            title = f'{tr[1].stats.network}.{tr[1].stats.station}.{tr[1].stats.location}.{tr[1].stats.channel} − starting {tr[1].stats["starttime"]}'                        
            torg = tr[1].times()
            if len(data) == 0:
                data = tr[0][:]
                fs = int(tr[0].stats.sampling_rate)
                title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
                torg = tr[0].times()
    except:
        try:
            p = "/scratch/irseppi/500sps/2019_0" + str(month) + "_" + str(day) + "/ZE_" + str(sta) + "_DPZ.msd"
            tr = obspy.read(p)

            tr.trim(tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs - wind, tr[0].stats.starttime + (int(h) * 3600) + (mins * 60) + secs + wind)

            data = tr[0][:]
            fs = int(tr[0].stats.sampling_rate)
            title = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'                        
            torg = tr[0].times()
        except:
            print(p)
            continue

    # Compute spectrogram
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
    # Error here with division by zero ##fix this
    spec, MDF = remove_median(Sxx)

    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    tprime0 = tarrive-start_time
    v0 = speed_mps
    height_m = alt - elev
    l = np.sqrt(dist_m**2 + (height_m)**2)

    tf = np.arange(0, 240, 1)

    coords = doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time,make_picks=False) 
    coords_array = np.array(coords)

    if len(coords) == 0:
        print('No picks for: ', date, flight_num, sta)
        continue
    # Convert the list of coordinates to a numpy array
    coords_array = np.array(coords)

    f0 = 116
    m0 = [f0, v0, l, tprime0, c]

    m,covm, F_m = invert_f(m0, coords_array, num_iterations=8)
    f0 = m[0]
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]
    
    ft = calc_ft(times, tprime0, f0, v0, l, c)
    if isinstance(sta, int):
        peaks = []
        p, _ = find_peaks(middle_column, distance = 7)
        corridor_width = (fs/2) / len(p) 
        if equip[0] == 'B' and equip[0:1] != 'BE':
            corridor_width = 3       
        if len(p) == 0:
            corridor_width = fs/4

        coord_inv = []

        for t_f in range(len(times)):
            upper = int(ft[t_f] + corridor_width)
            lower = int(ft[t_f] - corridor_width)
            if lower < 0:
                lower = 0
            if upper > len(frequencies):
                upper = len(frequencies)
            tt = spec[lower:upper, t_f]

            max_amplitude_index = np.argmax(tt)
            
            max_amplitude_frequency = frequencies[max_amplitude_index+lower]
            peaks.append(max_amplitude_frequency)
            coord_inv.append((times[t_f], max_amplitude_frequency))


        coord_inv_array = np.array(coord_inv)

        m,_,F_m = invert_f(m0, coord_inv_array,num_iterations=12)
        f0 = m[0]
        v0 = m[1]
        l = m[2]
        tprime0 = m[3]
        c = m[4]

        ft = calc_ft(times, tprime0, f0, v0, l, c)
        
        delf = np.array(ft) - np.array(peaks)
        
        new_coord_inv_array = []
        for i in range(len(delf)):
            if np.abs(delf[i]) <= 3:
                new_coord_inv_array.append(coord_inv_array[i])
        coord_inv_array = np.array(new_coord_inv_array)

        m,covm,F_m = invert_f(m0, coord_inv_array, num_iterations=12, sigma=5)
        
        f0 = m[0]
        v0 = m[1]
        l = m[2]
        tprime0 = m[3]
        c = m[4]

    mprior = []
    mprior.append(v0)
    mprior.append(l)
    mprior.append(tprime0) 
    mprior.append(c)      

    peaks, freqpeak =  overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight_num, sta, equip, closest_time, start_time, tprime0, wind, make_picks=False)

    f0_array = []
    w = len(peaks)
    for o in range(w):
        tprime = freqpeak[o]
        ft0p = peaks[o]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        mprior.append(f0)
        f0_array.append(f0)
    mprior = np.array(mprior)
                    
    corridor_width = 10
    if equip[0] == 'B' and equip[0:1] != 'BE':
        corridor_width = 3
    peaks_assos = []
    fobs = []
    tobs = []

    for pp in range(len(peaks)):
        tprime = freqpeak[pp]
        ft0p = peaks[pp]
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        ft = calc_ft(times,  tprime0, f0, v0, l, c)
        
        maxfreq = []
        coord_inv = []
        ttt = []

        f01 = f0 + corridor_width
        f02 = f0  - corridor_width
        upper = calc_ft(times,  tprime0, f01, v0, l, c)
        lower = calc_ft(times,  tprime0, f02, v0, l, c)

        for t_f in range(len(times)):

            try:      
                tt = spec[int(np.round(lower[t_f],0)):int(np.round(upper[t_f],0)), t_f]

                #For Boeing Jets
                if equip[0] == 'B' and equip[0:1] != 'BE':
                    max_amplitude_index,_ = find_peaks(tt, prominence = 5, wlen=5, height=vmax*0.5)
                else:
                    max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
                maxa = np.argmax(tt[max_amplitude_index])
                max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa])+int(np.round(lower[t_f],0))]
                maxfreq.append(max_amplitude_frequency)
                coord_inv.append((times[t_f], max_amplitude_frequency))
                ttt.append(times[t_f])

            except:
                continue

        if len(coord_inv) > 0:
            if f0 < 200:
                coord_inv_array = np.array(coord_inv)
                mtest = [f0,v0, l, tprime0,c]
                mtest,_, F_m = invert_f(mtest, coord_inv_array, num_iterations=4)
                ft = calc_ft(ttt,  mtest[3], mtest[0], mtest[1], mtest[2], c)
            else:
                ft = calc_ft(ttt,  tprime0, f0, v0, l, c)

            delf = np.array(ft) - np.array(maxfreq)

            count = 0
            for i in range(len(delf)):
                if np.abs(delf[i]) <= (4):
                    fobs.append(maxfreq[i])
                    tobs.append(ttt[i])
                    count += 1
            peaks_assos.append(count)


    if len(fobs) == 0:
        print('No picks for: ', date, flight_num, sta)
        continue

    try:
        tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, w, peaks_assos, make_picks=False)
        m, covm, f0_array, F_m = full_inversion(fobs, tobs, freqpeak, peaks, peaks_assos, tprime, tprime0, ft0p, v0, l, f0_array, mprior, c, w, 5)
    except:
        print('Inversion failed for: ', date, flight_num, sta)
        continue
    v0 = m[0]
    l = m[1]
    tprime0 = m[2]
    c = m[3]
    covm = np.sqrt(np.diag(covm))

    closest_index = np.argmin(np.abs(tprime0 - times))
    arrive_time = spec[:,closest_index]
    for i in range(len(arrive_time)):
        if arrive_time[i] < 0:
            arrive_time[i] = 0

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/with_c/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    qnum = plot_spectrgram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, covm, flight_num, middle_index, tarrive-start_time, closest_time, BASE_DIR, plot_show=False)

    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/with_c/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, BASE_DIR)
    
    #if rerun_fig == False:
    output.write(str(date)+','+str(flight_num)+','+str(sta)+','+str(closest_time)+','+str(tprime0)+','+str(v0)+','+str(l)+','+str(f0_array)+','+str(covm)+','+str(qnum)+','+str(Tc)+','+str(c)+','+str(F_m)+',\n') 

    #if rerun_fig == False:
    output.close()
