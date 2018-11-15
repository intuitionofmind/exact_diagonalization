#!/usr/bin/python3

import os
import math
import numpy as np
from scipy import linalg as la

from auxi import numSite, holeBasis, spinBasis, subDim, dim, numEval, numSam
from auxi import FindState, Flip

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

def LoadEigVal(f, paras):
    n = paras[2]
    rawData = np.fromfile(f % paras, dtype=np.float)
    eigVal = np.zeros((n, numEval))
    print(len(rawData))
    for i in range(n): 
        for j in range(numEval):
            eigVal[i][j] = rawData[i*numEval+j]
    return eigVal

def LoadEigVec(f, paras, l):
    rawData = np.fromfile(f % paras, dtype=np.float)
    eigVec = np.zeros((paras[0], dim), dtype=complex)
    for i in range(paras[0]):
        for j in range(dim):
            eigVec[i][j] = complex(rawData[l*(paras[0]*dim*2)+i*dim*2+j*2], rawData[l*(paras[0]*dim*2)+i*dim*2+j*2+1])
    return eigVec

def ChargeCurrent(v, hh, bond):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        s = spinBasis[i % subDim] #i = h*dim_sub+s
        h = holeBasis[i // subDim]
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'2'+b[h:] # To restore the t-J configuration. 
        if '0' != config[hh]:
            continue
        else:
            temp = list(config)
            k = bond[0]
            kk = bond[1]
            if h == bond[0]: # Hole's hopping from bond[0] to bond[1]. 
                hh = bond[1]
                temp[k], temp[kk] = temp[kk], temp[k]
                bb = ''
                for e in temp:
                    if e != '2':
                        bb = bb+e
                bb = bb[::-1]
                ss = int(bb, 2)
                l = FindState(ss, 0, subDim-1)
                sign = np.power(-1.0, k-kk+1)
                w[hh*subDim+l] = sign*v[i]
            elif h == bond[1]: # Hole's hopping from bond[1] to bond[0]. 
                hh = bond[0]
                temp[k], temp[kk] = temp[kk], temp[k]
                bb = ''
                for e in temp:
                    if e != '2':
                        bb = bb+e
                bb = bb[::-1]
                ss = int(bb, 2)
                l = FindState(ss, 0, subDim-1)
                sign = np.power(-1.0, kk-k+1)
                # w[hh*subDim+l] = 1.0*sign*v[i]  # For c_{i}^{\dagger}c_{j}+h.c.. 
                w[hh*subDim+l] = -1.0*sign*v[i] # For c_{i}^{\dagger}c_{j}-h.c.. 
    return w

dataDir = '/Users/wayne/Downloads/data/'
eigValsFile = 'eigenvalues%s_size1D%s_J0.3_10.0_%s_%s_sigma%s.dat'
eigVecsFile = 'eigenvectors%s_size1D%s_J0.3_10.0_%s_%s_sigma%s.dat'


# Set the slice in the array.
l = 0
J = 0.3

loadParas = (numEval, numSite, numSam, 'PBC', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), loadParas)
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile), loadParas, l)
print(loadParas, l, ene[l][:10])

sites = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
ps = 5
sites.remove(sp)
sitesLen = len(sites)
 
wf0 = wfArray[0]
wf1 = wfArray[1]
T = np.zeros((2, 2), dtype=complex)
if np.absolute(ene[l][0]-ene[l][1]) < 1e-10: # for degenerated eigenstates 
    wf0t = Translation(wf0)
    wf1t = Translation(wf1)
    T[0][0] = np.vdot(wf0, wf0t) # vdot() takes complex conjugate inner product
    T[0][1] = np.vdot(wf0, wf1t)
    T[1][0] = np.vdot(wf1, wf0t)
    T[1][1] = np.vdot(wf1, wf1t)
    w, v = la.eig(T) # w[i] is the eigenvalue of the corresponding vector v[:, i]
    print(w, np.angle(w, deg=True))
    wf00 = v[:, 0][0]*wf0+v[:, 0][1]*wf1
    wf11 = v[:, 1][0]*wf0+v[:, 1][1]*wf1
    if np.imag(w[0]) > 0:
        wf = wf00
    else:
        wf = wf11
else:
    wf = wf0
ccl = complex(0.0, 0.0)
bond = [0, 0]
for j in range(sitesLen):
    bond[0] = sites[j]
    bond[1] = sites[(j+1) % sitesLen]
    cc = np.vdot(wf, ChargeCurrentBond(wf, ps, bond))
    ccl += cc
    print(bond, format(cc), '.10f')
print(format(np.angle(w[0], deg=True), '.10f'), '  ', format(ccl, '.10f'))
