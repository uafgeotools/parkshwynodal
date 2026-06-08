from pyproj import Proj
from pathlib import Path
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from datetime import datetime, timezone

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

###############################################################

def make_base_dir(base_dir):

	"""
	Create a directory and its parent directories if they don't exist.

	Args:
		base_dir (str): The path of the directory to be created.

	Returns:
		None
	"""

	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path / parent
			if not current_path.exists():
				current_path.mkdir()

#########################################################################################################################################################################################################

def load_waveform(sta, arrive_time, spec_window=120, component="Z"):
	"""
	Load waveform data for a specific station and time window.

	Args:
		sta (str): Station code.
		arrive_time (float): Acoustic wave arrival time in seconds since the epoch.
		spec_window (int): Time window in seconds to trim the waveform data, default is 120 seconds.

	Returns:
		tuple: A tuple containing the waveform data, sampling frequency, time data correlating to waveform, and title
	"""
	
	ht = datetime.fromtimestamp(arrive_time, tz=timezone.utc)

	h = ht.hour
	mins = ht.minute
	secs = ht.second
	month = ht.month
	day = ht.day

	# Download waveform data from IRIS PH5WS
	client = Client(
		"http://service.iris.edu",
		service_mappings={
			"dataselect": "http://service.iris.edu/ph5ws/dataselect/1"
		}
	)
	starttime = UTCDateTime(2019, month, day, h, mins, secs) - spec_window
	endtime = starttime + spec_window*2

	tr = client.get_waveforms("ZE", str(sta), "*", f'DP{component}', starttime, endtime)
	tr = tr[0]
	data = tr.data
	sample_rate = int(tr.stats.sampling_rate)
	title = f'{tr.stats.network}.{tr.stats.station}.{tr.stats.location}.{tr.stats.channel} − starting {tr.stats["starttime"]}'
	t_wf = tr.times()


	return data, sample_rate, t_wf, title


#########################################################################################################################################################################################################