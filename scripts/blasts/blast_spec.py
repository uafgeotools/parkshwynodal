import obspy 
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram

tr = obspy.read("/scratch/naalexeev/NODAL/2019-03-02T01:00:00.000000Z.2019-03-02T02:00:00.000000Z.1225.mseed")

tr[2].trim(tr[2].stats.starttime + 42 * 60, tr[2].stats.starttime + 43 * 60 + 45)
data = tr[2][:]
fs = int(tr[2].stats.sampling_rate)
title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'                        
torg = tr[2].times()

#tr[2].spectrogram(per_lap=0.99,log=False,dbscale=True,cmap='hsv',clip=[0.4, 0.5])

# Compute spectrogram
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .99, detrend = 'constant') 

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     

ax1.plot(torg, data, 'k', linewidth=0.5)
ax1.set_title(title)

ax1.margins(x=0)
ax1.set_position([0.125, 0.6, 0.775, 0.3]) 

# Plot spectrogram
cax = ax2.pcolormesh(times, frequencies, Sxx, cmap='pink_r',  vmin = 0, vmax=7) #vmin=0.4*np.log(np.min(Sxx)), vmax=0.5*np.log(np.max(Sxx)))				
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Frequency (Hz)')

ax2.margins(x=0)
#ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

#plt.colorbar(mappable=cax, cax=ax3)
#ax3.set_ylabel('Relative Amplitude (dB)')

plt.show()