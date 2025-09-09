#after converting Alaska Railway kml file to text files reformat to be just lat and lon
#https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwi5tuT-_MP_AhUSJ30KHYmMDv0QFnoECAoQAQ&url=https%3A%2F%2Fwww.nps.gov%2Fdena%2Fplanyourvisit%2Fupload%2FAlaska%2520Railroad.kml&usg=AOvVaw2bKtne-33MDOx0zuXJjtR6

import sys
import fileinput

# replace all occurrences of ', 10' with '\n' 
for i, line in enumerate(fileinput.input('Alaska_Railroad.txt', inplace=1)):
    sys.stdout.write(line.replace(',', '   ')) 
    sys.stdout.write(line.replace(', 10 ', '\n')) 

