#!/usr/bin/python3

import os
import math
import cmath as cm
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from auxi import dim

def LoadEigVal(f, paras):
    n = paras[2]
    rawData = np.fromfile(f % paras, dtype=np.float)
    eigVal = np.zeros((n, numEval))
    # print(len(rawData))
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

# paras = (numEval, size, numSam, bouCon, 'False')
def TimeEvolutionUnitaryAppro(paras, sl, cutoff, t):
    ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
    wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile), paras, sl)
    TimeEvo = np.zeros((dim, dim), dtype=complex)
    for i in range(cutoff):
        Temp = np.outer(wfArray[i], np.conj(wfArray[i]))
        lam = ene[sl][i]-ene[sl][0]
        # print(lam)
        TimeEvo += np.exp(-1.j*lam*t)Temp
        # print(np.trace(TimeEvo))
        # print(i, TimeEvo[5][5], TimeEvo[5][100])
    TimeEvo = TimeEvo/np.sqrt(cutoff)
    return TimeEvo

J = 0.3

dataDir = '/Users/wayne/Downloads/data/'
eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_BC_%s_sigma%s.dat'
eigVecsFile = 'eigenvectors%s_size%s_J0.3_10.0_%s_BC_%s_sigma%s.dat'

bouCon = 'OO'
numEval = 279
size = 42
numSam = 5
s = 0
cutoff = 279 # The cut-off of eigenstate level to construct time evolution operator. 

paras = (numEval, size, numSam, bouCon, 'False')
np.set_printoptions(precision=5, suppress=True) # Only valid for printing numpy objects. 

wfArrayTest = LoadEigVec(os.path.join(dataDir, eigVecsFile), paras, 1)
wfInit = wfArrayTest[0]
for i in range(1):
    t = 0.01*i
    T = TimeEvolutionUnitaryAppro(paras, 0, cutoff, t)
    # print(np.vdot(T, T))
    # wfTime = np.dot(T, wfInit)
    # print(t, np.vdot(wfInit, wfTime))
