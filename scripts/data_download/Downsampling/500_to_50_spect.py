
import faulthandler
import numpy as np
import matplotlib.pyplot as plt
import obspy 

from obspy.signal.filter import bandpass
from tqdm import trange, tqdm

faulthandler.enable()

tr = obspy.read("/scratch/irseppi/nodal_data/500sps/2019_02_11/ZE_1001_DPZ.msd")
tr_new = obspy.read("/scratch/irseppi/nodal_data/50sps/2019_02_11/ZE_1001_DPZ.msd")

# Filtering with a bandpass on a copy of the original Trace
tr_filt = tr.copy()
tr_filt.filter('bandpass', freqmin=1, freqmax=25)

# Filtering with a bandpass on a copy of the downsampled Trace
tr_new_filt = tr_new.copy()
tr_new_filt.filter('bandpass', freqmin=1, freqmax=25)

for x in range(len(tr_filt)):

    t = np.arange(0, tr_filt[x].stats.npts / tr_filt[x].stats.sampling_rate, tr_filt[x].stats.delta)

    t_new = np.arange(0, tr_new_filt[x].stats.npts / tr_new_filt[x].stats.sampling_rate,tr_new_filt[x].stats.delta)
    
    plt.figure().set_figwidth(15)

    plt.plot(t, tr_filt[x].data, 'b', label='Raw', alpha=0.5)

    plt.plot(t_new[0:len(tr_new_filt[x])], tr_new_filt[x].data, 'r', label='Downsampled', alpha=0.4)

    plt.ylabel('Bandpassed Data')
    plt.xlabel('Time [s]')

    plt.suptitle(tr[x].stats.station+'.'+tr[x].stats.channel)
    plt.xlim(863,875)
    plt.legend()

    plt.show()
  

    t = np.arange(0, tr[x].stats.npts / tr[x].stats.sampling_rate, tr[x].stats.delta)

    t_new = np.arange(0, tr_new[x].stats.npts / tr_new[x].stats.sampling_rate,tr_new[x].stats.delta)

    plt.figure().set_figwidth(15)

    plt.plot(t, tr[x].data, 'b', label='Raw', alpha=0.5)

    plt.plot(t_new[0:len(tr_new[x])], tr_new[x].data, 'r', label='Downsampled', alpha=0.4)

    plt.ylabel('Bandpassed Data')
    plt.xlabel('Time [s]')

    plt.suptitle(tr[x].stats.station+'.'+tr[x].stats.channel)
    plt.xlim(863,875)
    plt.legend()

    plt.show()
  
