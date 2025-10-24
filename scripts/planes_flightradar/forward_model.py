import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy import sqrt as sympy_sqrt
import matplotlib as mpl

taken_derivative = False
if taken_derivative == True:
	f0, v0, tprime, tprime0, l, c = sp.symbols('f0 v0 tprime tprime0 l c')
	f = f0/(1+(v0/c)*(v0* ((tprime - tprime0 + l/c)- sympy_sqrt((tprime - tprime0 + l/c)**2-(1-v0**2/c**2)*((tprime - tprime0 + l/c)**2-l**2/c**2)))/(1-v0**2/c**2))/(sympy_sqrt(l**2+(v0* ((tprime - tprime0 + l/c)- sympy_sqrt((tprime - tprime0 + l/c)**2-(1-v0**2/c**2)*((tprime - tprime0 + l/c)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
	f = f0/(1+(v0/c)*(v0* ((tprime - tprime0)- sympy_sqrt((tprime - tprime0)**2-(1-v0**2/c**2)*((tprime - tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(sympy_sqrt(l**2+(v0* ((tprime - tprime0)- sympy_sqrt((tprime - tprime0)**2-(1-v0**2/c**2)*((tprime - tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
	# Take the derivative of f with respect to each variable
	for variable in [f0, v0, tprime0, l, c]:
		df = sp.diff(f, variable)
		print(f"Derivative of f with respect to {variable}:")
		sp.pprint(df)
		print("\n")

def get_t(v0,l,c,tprime0,tpr):
	t_array = []	
	for tprime in tpr:
		t = ((tprime - tprime0 + l/c)- np.sqrt((tprime - tprime0 + l/c)**2-(1-v0**2/c**2)*((tprime - tprime0 + l/c)**2-l**2/c**2)))/(1-v0**2/c**2)
		t_array.append(t)
	return t_array

def get_t_org(v0,l,c,tprime0,tpr):
	t_array = []	
	for tprime in tpr:
		t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
		t_array.append(t)
	return t_array

def get_f(v0,l,c,t_array,f0):
	ft = []
	for t in t_array:
		ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))	
		ft.append(ft0p)
	return ft

tpr = np.arange(0, 241, 1)

f0_vary = np.arange(0, 250, 20)
v0_vary = np.arange(0, 200, 20)
tprime0_vary = np.arange(0, 240, 20)
l_vary = np.arange(0, 5000, 200)
c_vary = np.arange(295, 355, 5)
c_vary = np.arange(200, 500, 5)
f0 = 150
tprime0 = 120
c = 320
v0 = 80
l = 4000
n = 1
test = False
if test:
	plt.figure(figsize=(10, 6))

	tprime = get_t(v0,l,c,tprime0,tpr)
	ft = get_f(v0,l,c,tprime, f0)

	plt.plot(tpr, ft, 'blue', linewidth=0.5, zorder=10)
	plt.axvline(tprime0, c='blue', linewidth=0.5, zorder=10)


	tprime = get_t_org(v0,l,c,tprime0,tpr)
	ft = get_f(v0,l,c,tprime, f0)

	plt.plot(tpr, ft, 'red', linewidth=0.5, zorder=10)
	plt.axvline(tprime0, c ='red', linewidth=0.5, zorder=10)

	plt.show()

fig, axs = plt.subplots(2, 3, figsize=(14, 10))
axs = axs.flatten()
for n in range(6):
	f0 = 150
	tprime0 = 120 
	c =  320
	v0 = 80
	l = 4000
	t0 = tprime0 - l/c
	tprime = get_t(v0, l, c, tprime0, tpr)
	ft = get_f(v0, l, c, tprime, f0)
	axs[n].plot(tpr, ft, 'k', linewidth=0.5, zorder=10)
	axs[n].axvline(tprime0, c='k', linewidth=0.5, zorder=10)
	if n == 0:
		axs[n].set_title('Base case')

	elif n == 1:
		axs[n].set_title('Varying v0')
		norm = plt.Normalize(np.min(v0_vary), np.max(v0_vary))
		cm = plt.cm.rainbow
		for v0_val in v0_vary:
			tprime = get_t(v0_val, l, c, tprime0, tpr)
			ft = get_f(v0_val, l, c, tprime, f0)
			axs[n].plot(tpr, ft, color=cm(norm(v0_val)), linewidth=0.5)
			t0 = tprime0 - l/c
			axs[n].axvline(t0, color=cm(norm(v0_val)), ls = '--', linewidth=0.5, zorder=10)
		sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
		sm.set_array([])
		cbar = plt.colorbar(sm, ax=axs[n], orientation='vertical', pad=0.02)

	elif n == 2:
		axs[n].set_title('Varying l')
		norm = plt.Normalize(np.min(l_vary), np.max(l_vary))
		cm = plt.cm.rainbow
		for l_val in l_vary:
			tprime = get_t(v0, l_val, c, tprime0, tpr)
			ft = get_f(v0, l_val, c, tprime, f0)
			axs[n].plot(tpr, ft, color=cm(norm(l_val)), linewidth=0.5)
			
		sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
		sm.set_array([])
		cbar = plt.colorbar(sm, ax=axs[n], orientation='vertical', pad=0.02)

	elif n == 3:
		axs[n].set_title('Varying c')
		norm = plt.Normalize(np.min(c_vary), np.max(c_vary))
		cm = plt.cm.rainbow
		for c_val in c_vary:
			tprime = get_t(v0, l, c_val, tprime0, tpr)
			ft = get_f(v0, l, c_val, tprime, f0)
			t0 = tprime0 - l/c_val
			axs[n].plot(tpr, ft, color=cm(norm(c_val)), linewidth=0.5)
		sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
		sm.set_array([])
		cbar = plt.colorbar(sm, ax=axs[n], orientation='vertical', pad=0.02)

	elif n == 4:
		axs[n].set_title("Varying t0'")
		norm = plt.Normalize(np.min(tprime0_vary), np.max(tprime0_vary))
		cm = plt.cm.rainbow
		for t0_val in tprime0_vary:
			tprime = get_t(v0, l, c, t0_val, tpr)
			ft = get_f(v0, l, c, tprime, f0)
			axs[n].plot(tpr, ft, color=cm(norm(t0_val)), linewidth=0.5)
		sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
		sm.set_array([])
		cbar = plt.colorbar(sm, ax=axs[n], orientation='vertical', pad=0.02)
		
	elif n == 5:
		axs[n].set_title('Varying f0')
		norm = plt.Normalize(np.min(f0_vary), np.max(f0_vary))
		cm = plt.cm.rainbow
		for f0_val in f0_vary:
			tprime = get_t(v0, l, c, tprime0, tpr)
			ft = get_f(v0, l, c, tprime, f0_val)
			axs[n].plot(tpr, ft, color=cm(norm(f0_val)), linewidth=0.5)
		sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
		sm.set_array([])
		cbar = plt.colorbar(sm, ax=axs[n], orientation='vertical', pad=0.02)
	axs[n].set_ylim(50, 250)
	axs[n].set_xlim(0, 240)
	if n == 0 or n == 3:
		axs[n].set_ylabel('Frequency (Hz)')
	if n == 3 or n == 4 or n == 5:
		axs[n].set_xlabel('Time (s)')

plt.tight_layout()
plt.show()
