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

def SpinCurrent(v, hh):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        # print(dim, i)
        s = spinBasis[i % subDim] #i = h*dim_sub+s
        h = holeBasis[i // subDim]
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        if h != hh:
            continue
        else:
            sites = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sites.remove(hh)
            ls = len(sites)
            for j in range(ls):
                k = sites[j]
                kk = sites[(j+1) % ls]
                if k < h and kk < h and b[k] != b[kk]:
                    ss = Flip(s, k, kk)
                    l = FindState(ss, 0, subDim-1)
                    sign = 1.0
                    if b[k] == '1':
                        sign = -1.0
                    w[h*subDim+l] += sign*v[i]
                elif k > h and kk < h and b[k-1] != b[kk]:
                    ss = Flip(s, k-1, kk)
                    l = FindState(ss, 0, subDim-1)
                    sign = 1.0
                    if b[k-1] == '1':
                        sign = -1.0
                    w[h*subDim+l] += sign*v[i]
                elif k < h and kk > h and b[k] != b[kk-1]:
                    ss = Flip(s, k, kk-1)
                    l = FindState(ss, 0, subDim-1)
                    sign = 1.0
                    if b[k] == '1':
                        sign = -1.0
                    w[h*subDim+l] += sign*v[i]
                elif k > h and kk > h and b[k-1] != b[kk-1]:
                    ss = Flip(s, k-1, kk-1)
                    l = FindState(ss, 0, subDim-1)
                    sign = 1.0
                    if b[k-1] == '1':
                        sign = -1.0
                    w[h*subDim+l] += sign*v[i]
    return w

dataDir = '/Users/wayne/Downloads/data/'
eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'
eigVecsFile = 'eigenvectors%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'

hp = 5

# Set the slice in the array.
l = 0
J = 0.3

loadParas = (numEval, numSite, numSam, 'PBC', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), loadParas)
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile), loadParas, l)
print(loadParas, l, ene[l][:40])

j = 0
for i in range(40):
    wf0 = wfArray[j]
    wf1 = wfArray[j+1]
    T = np.zeros((2, 2), dtype=complex)
    if np.absolute(ene[l][j]-ene[l][j+1]) < 1e-10: # for degenerated eigenstates 
        wf0t = Translation(wf0)
        wf1t = Translation(wf1)
        T[0][0] = np.vdot(wf0, wf0t) # vdot() takes complex conjugate inner product
        T[0][1] = np.vdot(wf0, wf1t)
        T[1][0] = np.vdot(wf1, wf0t)
        T[1][1] = np.vdot(wf1, wf1t)
        w, v = la.eig(T) # w[i] is the eigenvalue of the corresponding vector v[:, i]
        # print(j, np.angle(w, deg=True))
        wf00 = v[:, 0][0]*wf0+v[:, 0][1]*wf1
        wf11 = v[:, 1][0]*wf0+v[:, 1][1]*wf1
        sc0 = np.vdot(wf00, SpinCurrent(wf00, hp))
        sc1 = np.vdot(wf11, SpinCurrent(wf11, hp))
        print(j, '  ', np.around(ene[l][j], decimals=8), '  ', np.around(np.angle(w[0], deg=True), decimals=1), '  ', np.around(np.imag(sc0), decimals=8))
        print(j+1, '  ', np.around(ene[l][j+1], decimals=8), '  ', np.around(np.angle(w[1], deg=True), decimals=1), '  ', np.around(np.imag(sc1), decimals=8))
        j += 2
    else:
        wf = wf0
        sc = np.vdot(wf, SpinCurrent(wf, hp))
        print(j, '  ', np.around(ene[l][j], decimals=8), '  ', np.around(np.imag(sc), decimals=8))
        j += 1
 
# for sigma t-J model

loadParas = (numEval, numSite, numSam, 'PBC', 'True')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), loadParas)
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile), loadParas, l)
print(loadParas, l, ene[l][:40])

j = 0
for i in range(40):
    wf0 = wfArray[j]
    wf1 = wfArray[j+1]
    T = np.zeros((2, 2), dtype=complex)
    if np.absolute(ene[l][j]-ene[l][j+1]) < 1e-10: # for degenerated eigenstates 
        wf0t = Translation(wf0)
        wf1t = Translation(wf1)
        T[0][0] = np.vdot(wf0, wf0t) # vdot() takes complex conjugate inner product
        T[0][1] = np.vdot(wf0, wf1t)
        T[1][0] = np.vdot(wf1, wf0t)
        T[1][1] = np.vdot(wf1, wf1t)
        w, v = la.eig(T) # w[i] is the eigenvalue of the corresponding vector v[:, i]
        # print(j, np.angle(w, deg=True))
        wf00 = v[:, 0][0]*wf0+v[:, 0][1]*wf1
        wf11 = v[:, 1][0]*wf0+v[:, 1][1]*wf1
        sc0 = np.vdot(wf00, SpinCurrent(wf00, hp))
        sc1 = np.vdot(wf11, SpinCurrent(wf11, hp))
        print(j, '  ', np.around(ene[l][j], decimals=8), '  ', np.around(np.angle(w[0], deg=True), decimals=1), '  ', np.around(np.imag(sc0), decimals=8))
        print(j+1, '  ', np.around(ene[l][j+1], decimals=8), '  ', np.around(np.angle(w[1], deg=True), decimals=1), '  ', np.around(np.imag(sc1), decimals=8))
        j += 2
    else:
        wf = wf0
        sc = np.vdot(wf, SpinCurrent(wf, hp))
        print(j, '  ', np.around(ene[l][j], decimals=8), '  ', np.around(np.imag(sc), decimals=8))
        j += 1
