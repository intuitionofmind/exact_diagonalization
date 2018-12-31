import numpy as np

numSite = 8

e = 0.0
for m in range(-3, 3, 1):
    print(m)
    k = 2*np.pi*m / numSite
    e += -2.0*np.cos(k)

print(e / numSite)
