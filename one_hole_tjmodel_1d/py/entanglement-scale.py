#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes

fe = '/home/zhengwei/Academics/physics/projects/exact_diag_tj_model/data/totSZ_half/eigenvalues_%s_len_%s_J-range_0.0-0.1-%s_%s_sigma_%s_flux_%s-PI.dat'
fv = '/home/zhengwei/Academics/physics/projects/exact_diag_tj_model/data/totSZ_half/eigenvectors_%s_len_%s_J-range_0.0-0.1-%s_%s_sigma_%s_flux_%s-PI.dat'

fs = 16

def FindState(s, lower, upper):
    if lower == upper:
        return upper
    else:
        mid = (lower+upper)//2
        if s > spinBasis[mid]:
            return FindState(s, mid+1, upper)
        else:
            return FindState(s, lower, mid)

def Flip(s, i, j):
    f = int(np.power(2, i))+int(np.power(2, j))
    return s^f

def Translation(v):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        h = holeBasis[i // subDim]
        hh = (h+1) % numSite
        if h == numSite-1:
            w[hh*subDim+(i % subDim)] = v[i]
        else:
            s = spinBasis[i % subDim]
            b = bin(s)[2:]
            b = b.rjust(numSite, '0')
            b = b[::-1]
            bb = b[numSite-2]+b[:numSite-2]
            bb = bb[::-1]
            ss = int(bb, 2)
            l = FindState(ss, 0, subDim-1)
            w[hh*subDim+l] = np.power(-1.0, numSite)*v[i]
    return w
 
def LoadEigVal(paras):
    rawData = np.fromfile(fe % paras, dtype=np.float)
    eigVal = np.zeros((paras[2], paras[0]))
    for i in range(paras[2]): 
        for j in range(paras[0]):
            eigVal[i][j] = rawData[i*paras[0]+j]
    return eigVal

def LoadEigVec(paras, k):
    rawData = np.fromfile(fv % paras, dtype=np.float)
    eigVec = np.zeros((paras[2], dim), dtype=complex)
    for i in range(paras[2]): 
        for j in range(dim):
            eigVec[i][j] = complex(rawData[i*(paras[0]*dim*2)+k*dim*2+j*2], rawData[i*(paras[0]*dim*2)+k*dim*2+j*2+1])
    return eigVec

def EntanglementSpectrum(wf):
    Rho = np.zeros((paras[1], paras[1]), dtype=complex)
    for i in range(0, paras[1]):
        for j in range(0, paras[1]):
            for k in range(0, subDim):
                Rho[i][j] += wf[i*subDim+k]*np.conjugate(wf[j*subDim+k])
    np.set_printoptions(precision=3, suppress=True)
    es = np.linalg.eigvalsh(Rho)
    return es

def EntanglementEntropy(wf):
    es = EntanglementSpectrum(wf)
    ee = -1.0*np.dot(es, np.log(np.absolute(es)))  # entanglement entropy
    return ee

m = 8 # m samples along J/t 
n = 7 # for a specific J/t, there are n sizes  
numEval = 5
numSam = 10
flux = 0.0
mz = 0.5

x = np.linspace(8, 8+(n-1)*2, num=n)
xx = np.log(x)
y = np.zeros(n)

'''
EE = np.zeros((m, n))

for i in range(n):
    numSite = 8+i*2
    paras = (numEval, numSite, numSam, 'PBC', 'false', flux)

    # Hilbert space construction
    max = int(math.pow(2, numSite-1))
    n = 0
    spinBasis = np.zeros(0, dtype=int)
    holeBasis = np.zeros(numSite, dtype=int)
    for ii in range(max):
        if bin(ii).count('1') == int(mz+(numSite-1)/2.0):
            spinBasis = np.append(spinBasis, ii)
            n += 1
    for ii in range(numSite):
        holeBasis[ii] = ii
    subDim = n
    dim = numSite*n

    ene = LoadEigVal(paras)
    wf0 = LoadEigVec(paras, 0)
    wf1 = LoadEigVec(paras, 1)
    for j in range(m):
        T = np.zeros((2, 2), dtype=complex)
        if np.absolute(ene[j][0]-ene[j][1]) < 1e-10:
            wf0t = Translation(wf0[j])
            wf1t = Translation(wf1[j])
            T[0][0] = np.vdot(wf0[j], wf0t)
            T[0][1] = np.vdot(wf0[j], wf1t)
            T[1][0] = np.vdot(wf1[j], wf0t)
            T[1][1] = np.vdot(wf1[j], wf1t)
            w, v = la.eig(T)
            if w[0].imag > 0:
                choice = 0
            else:
                choice = 1
            # print(i, np.angle(w, deg=True))
            wf = v[:, choice][0]*wf0[j]+v[:, choice][1]*wf1[j] 
        else:
            wf = wf0[j]
        EE[j][i] = EntanglementEntropy(wf)
        print(numSite, j)

f = 'ee-scale'
np.save(f, EE)
'''
EE = np.load('ee-scale.npy')
print(EE)

fig = plt.figure(figsize=(8, 6))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax1 = plt.subplot(111)
p1 = np.polyfit(xx, EE[1], 1)
p2 = np.polyfit(xx, EE[3], 1)
p3 = np.polyfit(xx, EE[5], 1)
p4 = np.polyfit(xx, EE[7], 1)
print(p1[0], p2[0], p3[0], p4[0])
# plt.plot(xx, EE[1], 'o', color='red', label='$J/t=0.1$')
# plt.plot(xx, xx*p1[0]+p1[1], '--', color='red')
# plt.plot(xx, EE[3], '>', color='blue', label='$J/t=0.3$')
# plt.plot(xx, xx*p2[0]+p2[1], '--', color='blue')
# plt.plot(xx, EE[5], 'd', color='green', label='$J/t=0.5$')
# plt.plot(xx, xx*p3[0]+p3[1], '--', color='green')
# plt.plot(xx, EE[7], 'x', color='cyan', label='$J/t=0.7$')
# plt.plot(xx, xx*p4[0]+p4[1], '--', color='cyan')

plt.plot(xx, EE[2], 'o', color='red', linestyle='--', label='$J/t=0.2$')
plt.plot(xx, EE[4], '>', color='blue', linestyle='--', label='$J/t=0.4$')
plt.plot(xx, EE[6], 'd', color='green', linestyle='--', label='$J/t=0.6$')
# plt.plot(xx, EE[7], 'x', color='cyan', linestyle='--', label='$J/t=0.7$')


plt.xlabel('$N$', fontsize=fs)
plt.ylabel('Mutual spin-charge EE', fontsize=fs)
plt.xlim([2.07, 3.0])
# plt.xlim([-0.15, 0.5])
plt.legend(loc='best', frameon=False, fontsize=fs)
plt.xticks(np.log([8, 10, 12, 14, 16, 18, 20]), [8, 10, 12, 14, 16, 18, 20], fontsize=fs)
plt.setp(ax1.get_yticklabels(), fontsize=fs)
ax1.text(0.4, 0.8, '$\sigma\cdot{t}$-$J$', fontsize=fs, transform=ax1.transAxes)

image = 'ee-scale.png'
fig.tight_layout()
plt.savefig(image, format='PNG')
# plt.show()
