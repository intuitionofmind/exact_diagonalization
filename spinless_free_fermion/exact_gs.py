import numpy as np

numSite = 16
numEle = 8

e = 0.0
for m in range(-3, 5, 1):
    print(m)
    k = 2*np.pi*m / numSite
    e += -2.0*np.cos(k)

print(e / numSite)
