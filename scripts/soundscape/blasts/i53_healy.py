from waveform_collection import gather_waveforms
from obspy.core import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from array_processing.tools.plotting import array_plot
from lts_array import ltsva

# Data collection parameters
SOURCE = 'IRIS'
NETWORK = 'IM'
STATION = 'I53H?'
LOCATION = '*'
CHANNEL = 'BDF'
START = UTCDateTime('2019-04-13 00:19')
END = UTCDateTime('2019-04-13 00:29')

# Filtering
FMIN = 1  # [Hz]
FMAX = 3  # [Hz]

# Array processing
WINLEN = 30  # [s]
WINOVER = 0.9

# Grab and filter waveforms
st = gather_waveforms(
    SOURCE, NETWORK, STATION, LOCATION, CHANNEL, START, END, remove_response=True
)
st.filter('bandpass', freqmin=FMIN, freqmax=FMAX, corners=2, zerophase=True)
st.taper(max_percentage=0.01)

latlist = [tr.stats.latitude for tr in st]
lonlist = [tr.stats.longitude for tr in st]

# Array process
vel, baz, t, mdccm, stdict, sigma_tau, conf_int_vel, conf_int_baz = ltsva(
    st, latlist, lonlist, WINLEN, WINOVER, alpha=1
)

# Plot
fig, axs = array_plot(st, t, mdccm, vel, baz, ccmplot=True)
true_baz = gps2dist_azimuth(latlist[0], lonlist[0], 63.99, -148.74)[1]
axs[3].axhline(true_baz, zorder=-5, color='gray', linestyle='--')
fig.show()
