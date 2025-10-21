import faulthandler

from obspy import read
from pathlib import Path

faulthandler.enable()

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()


for day in range (14,22):
	BASE_DIR = "/home/irseppi/nodal_data/50sps/2019_02_"+str(day)
	make_base_dir(BASE_DIR)

	for z in range(0,3):
		for y in range(9,10):
			try:
				tr = read("/home/irseppi/nodal_data/500sps/2019_02_"+str(day)+"/ZE_1"+str(z)+str(y)+"*.msd")
				print('READ')
				
				tr_new = tr.copy()
				tr_new.resample(50)

				for x in range(len(tr_new)):
					save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
					print(save_name)

					save_path = (Path(BASE_DIR) / save_name)

					tr_new[x].write(str(save_path), format='MSEED')
					print('SAVED')
			except:
				continue
	try:
			tr = read("/home/irseppi/nodal_data/500sps/2019_02_"+str(day)+"/ZE_130*.msd")
			print('READ')
			
			tr_new = tr.copy()
			tr_new.resample(50)

			for x in range(len(tr_new)):
				save_name = BASE_DIR + "/ZE_{}_{}.msd".format(tr_new[x].stats.station,tr_new[x].stats.channel)
				print(save_name)

				save_path = (Path(BASE_DIR) / save_name)

				tr_new[x].write(str(save_path), format='MSEED')
				print('SAVED')
	except:
		continue


