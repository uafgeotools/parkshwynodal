import sys
import fileinput
import matplotlib.pyplot as plt
import numpy as np

# explicit function to return the letter count
def numFrequency(fileName, num):
	# open file in read mode
	text = open(fileName, "r")

	# declare count variable
	count = 0
	for line in text.readlines():
		val = line.split()
		# compare each character with
		if int(val[0]) == int(num):
			count += 1
		else:
			continue

	# return letter count
	return count

lab=[]	
num=[]
for x in range(1001, 1306):
	# call function and display the letter count
	count=numFrequency('nodalt.arrival', x)
	if count != 0:
		lab.append(str(x))
		num.append(int(count))
	else:
		continue

fig, ax = plt.subplots()

ax.bar(lab, num)

ax.set_ylabel('Number of Arrivals')
ax.set_title('Number of Arrival Picks by Station')
ax.set_xticklabels(lab,rotation=45)

plt.show()

#if count text file is already dowloaded:
text = open('count.txt', "r")
lab=[]	
num=[]
for line in text.readlines():
	val = line.split()
	lab.append(str(val[0]))
	num.append(int(val[1]))

fig, ax = plt.subplots()

ax.bar(lab, num)

ax.set_ylabel('Number of Arrivals')
ax.set_title('Number of Arrival Picks by Station')
ax.set_xticklabels(lab,rotation=90)

plt.show()

