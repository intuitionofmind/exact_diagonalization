#!/usr/bin/python3

import os
import math
import cmath as cm
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from auxi import holeBasis, spinBasis, subDim, dim, numSiteX, numSiteY, numEval, numSam
from auxi import FindState, Flip

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

def ChargeCurrentBond(v, proSp, bond):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        s = spinBasis[i % subDim] # i = h*dim_sub+s
        h = holeBasis[i // subDim]
        if h not in bond:
            continue
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'2'+b[h:] # To restore the t-J configuration. 
        if '0' != config[proSp]:
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
                # print(h, b, hh, bb)
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
                w[hh*subDim+l] = 1.0*sign*v[i]  # For c_{i}^{\dagger}c_{j}+h.c.. 
                # w[hh*subDim+l] = -1.0*sign*v[i] # For c_{i}^{\dagger}c_{j}-h.c.. 
    return w

dataDir = '/Users/wayne/Downloads/data/'
eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_BC_%s_sigma%s.dat'
eigVecsFile = 'eigenvectors%s_size%s_J0.3_10.0_%s_BC_%s_sigma%s.dat'

numSite = numSiteX*numSiteY

size = 44
l = 0
J = 0.3
fold = 6

paras = (numEval, size, numSam, 'PP', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
print(paras, l, ene[l][:fold])

f = './translationWaveFunction_fold%s_J%s_%s_sigma%s.npy'
npyParas = (fold, J, 'PP', 'False')
wfArray = np.load(f % npyParas)
# wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile), paras, l)

proSpin = 5
sites = [0, 4, 8, 9, 10, 6, 2, 1]
sitesLen = len(sites)
for i in range(fold):
    wf = wfArray[i]
    ccl = complex(0.0, 0.0)
    bond = [0, 0]
    for j in range(sitesLen):
        bond[0] = sites[j]
        bond[1] = sites[(j+1) % sitesLen]
        cc = np.vdot(wf, ChargeCurrentBond(wf, proSpin, bond))
        print(bond, '  ', format(cc, '.10f'))
        ccl += cc
    print('summation along the loop', format(ccl, '.10f'))
