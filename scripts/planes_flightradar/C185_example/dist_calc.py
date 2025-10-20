import numpy as np

tprime = 23.92
t0 = 116.18
v0 = 69.45
l = 952.18
c = 313

t = ((tprime - t0)- np.sqrt((tprime - t0)**2-(1-v0**2/c**2)*((tprime - t0)**2-l**2/c**2)))/(1-v0**2/c**2)
print(v0*t)
