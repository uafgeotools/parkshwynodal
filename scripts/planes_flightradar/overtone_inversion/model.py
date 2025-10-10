import numpy as np
import matplotlib.pyplot as plt

tpr = np.arange(0, 241, 1)

f0_vary = np.arange(0, 250, 20)
v0_vary = np.arange(0, 200, 20)
t0_vary = np.arange(0, 250, 24)
l_vary = np.arange(0, 2000, 250)
c_vary = np.arange(295, 355, 5)
for n in range(3, 4):
	f0 = 100
	t0 = 120
	c = 343
	v0 = 50
	l = 1000

	plt.figure(figsize=(10, 6))
	ft = []
	if n == 0:
		fchange = []
		#norm = plt.Normalize(np.min(f0_vary), np.max(f0_vary))
		#cm = plt.cm.rainbow

		#sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
		plt.title('Varying f0')
		for f0 in f0_vary:

			ft = []
			for tprime in tpr:
				ft0p = f0*1/(1+(v0/c)*(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
				
				ft.append(ft0p)
			fchange.append(ft[-1]-ft[0])
			#plt.plot(tpr, ft, color=cm(norm(v0)),  linewidth=0.5)
		plt.plot(f0_vary, fchange) #,  color=cm(norm(f0)), linewidth=0.5)
		plt.xlabel('f0')

	if n == 1:
		fchange = []
		#norm = plt.Normalize(np.min(v0_vary), np.max(v0_vary))
		#cm = plt.cm.rainbow

		#sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
		plt.title('Varying v0')
		for v0 in v0_vary:
			ft = []
			for tprime in tpr:
				ft0p = f0*1/(1+(v0/c)*(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
			
				ft.append(ft0p)
			fchange.append(ft[-1]-ft[0])
			#plt.plot(tpr, ft, color=cm(norm(v0)),  linewidth=0.5)
		plt.plot(v0_vary, fchange)
		plt.xlabel('v0')
	if n == 2:
		#norm = plt.Normalize(np.min(l_vary), np.max(l_vary))
		#cm = plt.cm.rainbow

		#sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
		plt.title('Varying l')
		
		fchange = []
		for l in l_vary:
			ft = []
			for tprime in tpr:
				ft0p = f0*1/(1+(v0/c)*(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
				
				ft.append(ft0p)
			
			fchange.append(ft[-1]-ft[0])
			#plt.plot(tpr, ft, color=cm(norm(l)), linewidth=0.5)
		plt.plot(l_vary, fchange)
		plt.xlabel('l')
	if n == 3:
		fchange = []
		norm = plt.Normalize(np.min(c_vary), np.max(c_vary))
		cm = plt.cm.rainbow

		sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
		plt.title('Varying c')
		for c in c_vary:
			ft = []
			for tprime in tpr:
				ft0p = f0*1/(1+(v0/c)*(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
			
				ft.append(ft0p)
			fchange.append(ft[-1]-ft[0])
			plt.plot(tpr, ft, color=cm(norm(c)),  linewidth=0.5)

		#plt.plot(c_vary, fchange)
		#plt.xlabel('c')
	plt.ylim(83,123)
	plt.ylabel('Frequency (Hz)')
	#plt.show()

#plt.figure()
tpr = np.arange(0, 241, 1)
f0 = 100
t0 = 120
c = 343
v0 = 50
l = 1000
ft = []

for tprime in tpr:
	ft0p = f0*1/(1+(v0/c)*(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
			
	ft.append(ft0p)
plt.plot(tpr, ft, 'k', linewidth=0.5)
l = 10
ti = [50,200]
for t0 in ti:
	
	ft0p = f0*1/(1+(v0/c)*(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - t0)- np.sqrt((tprime-t0)**2-(1-v0**2/c**2)*((tprime-t0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
	#plt.axhline(ft0p, color='k', linestyle='--')

plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.xlim(0, 240)
#plt.ylim(0, 250)
plt.show()					
