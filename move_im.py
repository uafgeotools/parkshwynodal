import os
import shutil
#open files in 
BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/sigma/'
BASE_DIR2 = '/scratch/irseppi/nodal_data/plane_info/sigma_comp/'

for file_name in os.listdir(BASE_DIR):
    sig_file = os.path.join(BASE_DIR, file_name)
    sig_hold = os.path.join(BASE_DIR2, file_name)
    for eq in os.listdir(sig_file):
        eq_file = os.path.join(sig_file, eq)
        for date_file in os.listdir(eq_file):
            full_path = os.path.join(eq_file, date_file)
            for flight_file in os.listdir(full_path):
                flight_path = os.path.join(full_path, flight_file)
                for sta_file in os.listdir(flight_path):
                    sta_path = os.path.join(flight_path, sta_file)
                    for file in os.listdir(sta_path):
                        image_file = os.path.join(sta_path, file)
                        new_image_file = os.path.join(sig_hold,file)
                        os.makedirs(os.path.dirname(new_image_file), exist_ok=True)
                        shutil.copy2(image_file, new_image_file)
