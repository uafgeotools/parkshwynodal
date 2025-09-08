import gc

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from pathlib import Path
from threading import Thread
import math


# Writen by Nicholas Alexeev 
# From https://github.com/scifi6546/denali_nodal_set/blob/master/par_download.py 

def make_base_dir(base_dir):
    base_dir = Path(base_dir)
    if not base_dir.exists():
        current_path = Path("/")
        for parent in base_dir.parts:
            current_path = current_path / parent
            if not current_path.exists():
                current_path.mkdir()


def thread_task(id):
    time_difference = END_DATE - START_DATE
    thread_time_window = time_difference / float(NUM_THREADS)
    start_offset = thread_time_window * id
    end_offset = start_offset + thread_time_window
    start_time = START_DATE + start_offset
    end_time = START_DATE + end_offset
    print_str = "thread id: {}\n\tthread download start: {} end: {}\n".format(id, start_time, end_time)
    print(print_str)
    download_waveforms(start_time, end_time, DOWNLOAD_WINDOW,id)


def download_waveforms(download_start_time, download_end_time, download_window,thread_id):
    waveform_client = Client("http://service.iris.edu",
                             service_mappings={"dataselect": "http://service.iris.edu/ph5ws/dataselect/1"},
                             user="ctape@alaska.edu", password="dln3mjKtap3m9")
    diff = download_end_time - download_start_time
    for i in range(0, int(math.ceil(diff / download_window))):
        download_start = download_start_time + float(i) * download_window
        download_end = download_end_time + float(i + 1) * download_window
        for station in stations_to_download:
            print_message = ""

            save_name = BASE_DIR + "{}.{}.{}.mseed".format(download_start, download_end, station)
            print("(thread: {}) downloading file: {}".format(thread_id,save_name))
            print_message += "downloading data from station {}, {} to {}\n".format(station, download_start,
                                                                                   download_end)

            try:
                if not Path(save_name).exists():
                    waveform = waveform_client.get_waveforms(network="ZE", location="*", station=station, channel="*",
                                                             starttime=download_start,
                                                             endtime=download_end)
                    print_message += "\tsaving to file {}".format(save_name)

                    waveform.write(save_name)
                else:
                    print_message += "file {} already exists, skipping".format(save_name)

            except Exception as e:
                print_message += "download error for file {}: {}".format(save_name, e)
            print(print_message)

            gc.collect()


BASE_DIR = "/scratch/naalexeev/NODAL/"
NUM_THREADS = 1
make_base_dir(BASE_DIR)

START_DATE = UTCDateTime(2019, 2, 13)
END_DATE = UTCDateTime(2019, 3, 26)
# diff = end_date-start_date
DOWNLOAD_WINDOW = 60.0 * 60.0
STATION = 'http://service.iris.edu/ph5ws/station/1'
c = Client("http://service.iris.edu", service_mappings={'station': STATION}, debug=False)
stations = c.get_stations(network="ZE", location="*", station="*", channel="*",
                          starttime=UTCDateTime("2019-02-25T18:20:50.906000Z"),
                          endtime=UTCDateTime("2019-02-25T18:27:30.906000Z"), minlatitude=None, maxlatitude=None,
                          minlongitude=None, maxlongitude=None, level="response")
stations.write(BASE_DIR + "/stations.xml", format="STATIONXML")

# rint(stations)
stations_to_download = []
for net in stations:
    for s in net:
        stations_to_download.append(s.code)

print("downloading data from {} to {} for stations: ".format(START_DATE, END_DATE))
for station in stations_to_download:
    print("\t{}".format(station))

for i in range(0, NUM_THREADS):
    thread_task(i)

