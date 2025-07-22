import os
def delete_empty_file(file_path):
    """
    Checks if a given file is empty and deletes it if it is.

    Args:
        file_path (str): The path to the file to be checked.
    """
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' does not exist.")
        return

    try:
        if os.path.getsize(file_path) == 0:
            os.remove(file_path)
            print(f"File '{file_path}' was empty and has been deleted.")
        else:
            print(f"File '{file_path}' is not empty.")
    except OSError as e:
        print(f"Error processing file '{file_path}': {e}")
        
# loop through the directories in the directory
for dir_name in os.listdir('/home/irseppi/REPOSITORIES/parkshwynodal/output/'):
    dir_path = os.path.join('/home/irseppi/REPOSITORIES/parkshwynodal/output/', dir_name)
    if  os.path.isdir(dir_path):
        for equip_dir in os.listdir(dir_path):
            equip_dir_path = os.path.join(dir_path, equip_dir)
            if os.path.isdir(equip_dir_path):
                for date_dir in os.listdir(equip_dir_path):
                    date_dir_path = os.path.join(equip_dir_path, date_dir)
                    if os.path.isdir(date_dir_path):
                        for flight_dir in os.listdir(date_dir_path):
                            flight_dir_path = os.path.join(date_dir_path, flight_dir)
                            if os.path.isdir(flight_dir_path):
                                for sta_dir in os.listdir(flight_dir_path):
                                    sta_dir_path = os.path.join(flight_dir_path, sta_dir)
                                    if os.path.isdir(sta_dir_path):
                                        for file_name in os.listdir(sta_dir_path):
                                            file_path = os.path.join(sta_dir_path, file_name)
                                            delete_empty_file(file_path)

