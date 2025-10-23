from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from pathlib import Path
import sqlite3
import math

# Writen by Nicholas Alexeev 
# From https://github.com/scifi6546/denali_nodal_set/blob/master/download_gaps.py 

def make_base_dir(base_dir):
    base_dir = Path(base_dir)
    if not base_dir.exists():
        current_path = Path("/")
        for parent in base_dir.parts:
            current_path = current_path/parent
            if not current_path.exists():
                current_path.mkdir()



BASE_DIR = "/scratch/naalexeev/NODAL/"
make_base_dir(BASE_DIR)

old_start_date =UTCDateTime(2019,2,13)
end_date = UTCDateTime(2019,3,26)
start_date = old_start_date + (end_date-old_start_date)/2.0

diff = end_date-start_date
DOWNLOAD_WINDOW = 60.0 * 60.0


conn = sqlite3.connect("/scratch/naalexeev/flight_database.sqlite")
num_gaps = conn.execute("SELECT COUNT(station_id) FROM waveforms WHERE path IS NULL;").fetchall()
print("num gaps: {}".format(num_gaps))
gaps_rows = conn.execute("SELECT DISTINCT station_id,start_time,end_time FROM waveforms WHERE path IS NULL;").fetchall()
gaps_data = []
for gap in gaps_rows:
    gaps_data.append({"station_id":gap[0],"starttime":UTCDateTime(gap[1]),"endtime":UTCDateTime(gap[2]),"start_timestamp":gap[1],"end_timestamp":gap[2]})
waveform_client = Client("http://service.iris.edu",
                   service_mappings={"dataselect": "http://service.iris.edu/ph5ws/dataselect/1"})
i=0
for gap in gaps_data:
    save_name = BASE_DIR + "{}.{}.{}.mseed".format(gap["starttime"], gap["endtime"], gap["station_id"])
    save_path = (Path(BASE_DIR) / save_name)
    save_str = str(save_path)
    sql_str = "UPDATE waveforms SET path = '{}' WHERE station_id = '{}' AND start_time = {} and end_time = {}".format(
        str(save_str),
        gap["station_id"],
        gap["start_timestamp"],
        gap["end_timestamp"]
    )



    if not save_path.exists():
        print("downloading file:{}\n\t{}%\n\tsql exec str: {}".format(str(save_path),100.0*float(i)/float(len(gaps_data)),sql_str))


        try:
            waveform = waveform_client.get_waveforms(network="ZE", location="*", station=gap["station_id"], channel="*",
                                          starttime=gap["starttime"],
                                          endtime=gap["endtime"])
            waveform.write(str(save_path))
            conn.execute(sql_str)
            conn.commit()
        except Exception as e:
            print(e)
    else:
        conn.execute(sql_str)
        conn.commit()
    i+=1

