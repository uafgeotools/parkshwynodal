#after converting 3 Parks Highway kml(https://dot.alaska.gov/stwddes/gis/transdata/GE_Files/Alaska_DOTPF_RoadSystem_081114.kml) files to text files reformat to be just lat and lon
import sys
import fileinput

for i, line in enumerate(fileinput.input('PARKSHIGHWAY3.txt', inplace=1)):
    sys.stdout.write(line.replace(',0', '\n')) 
    sys.stdout.write(line.replace('</coordinates>..<coordinates>', '  '))

