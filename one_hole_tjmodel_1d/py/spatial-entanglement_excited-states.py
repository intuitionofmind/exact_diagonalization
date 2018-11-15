#!/usr/bin/python3

import os
import math
import numpy as np
from scipy import linalg as la
from auxi import *

cut = 5

# Count the dimension of the Hilbert space for subsystems.
basisA = []
basisB = []
indexDict = [] # Construct the dictonary which maps indices between two representations. 
for i in range(numSite):
    h = holeBasis[i]
    for j in range(subDim):
        s = spinBasis[j]
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'2'+b[h:] # To restore the t-J configuration. 
        configA = config[:cut] 
        configB = config[cut:]
        # print(config, configA, configB)
        if configA not in basisA:
            basisA.append(configA)
        if configB not in basisB:
            basisB.append(configB)
        indexDict.append((basisA.index(configA), basisB.index(configB)))

dimA = len(basisA)
dimB = len(basisB)

print(dim)
print(dimA, dimB)

def Reshape(wf):
    WF = np.zeros((dimA, dimB), dtype=complex)
    for i in range(dim):
        '''
        s = spinBasis[i % subDim] # i = h*sumDim+(i % subDim)
        h = holeBasis[i // subDim]
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'2'+b[h:] # To restore the t-J configuration. 
        configA = config[:cut] 
        configB = config[cut:]
        j = basisA.index(configA)
        k = basisB.index(configB)
        '''
        j = indexDict[i][0]
        k = indexDict[i][1]
        WF[j][k] = wf[i]
    return WF

def RedDenMatrixA(WF):
    Rho = np.zeros((dimA, dimA), dtype=complex)
    for i in range(dimA):
        for j in range(dimA):
            for k in range(dimB):
                Rho[i][j] += WF[i][k]*np.conjugate(WF[j][k])
    return Rho

def SpatialES(wf):
    WF = Reshape(wf)
    Rho = RedDenMatrixA(WF)
    es = np.linalg.eigvalsh(Rho)
    # print(es)
    return es

def SpatialEE(wf):
    es = SpatialES(wf)
    ee = -1.0*np.dot(es, np.log(np.absolute(es)))  # entanglement entropy
    print(ee)
    return ee

def SpatialEEArray(wfArray):
    eeArray = np.zeros(numEval)
    for i in range(numEval):
        eeArray[i] = SpatialEE(wfArray[i])
        print(i)
    return eeArray

def SpatialESArray(wfArray):
    esArray = np.zeros((numEval, dimA))
    for i in range(numEval):
        esArray[i] = SpatialES(wfArray[i])
    return esArray 

l = 9
fe = 'spatial-ee_len-%s_J-%s_%s_sigma-%s'
fs = 'spatial-es_len-%s_J-%s_%s_sigma-%s'

loadParas = (numEval, numSite, numSam, 'PBC', 'false')
saveParas = (numSite, (l+1)*step, 'PBC', 'false')
wfArray = LoadEigVec2(loadParas, l)
ee = SpatialEEArray(wfArray)
es = SpatialESArray(wfArray)
np.save(fe % saveParas, ee)
np.save(fs % saveParas, es)

loadParas = (numEval, numSite, numSam, 'PBC', 'true')
saveParas = (numSite, (l+1)*step, 'PBC', 'true')
wfArray = LoadEigVec2(loadParas, l)
ee = SpatialEEArray(wfArray)
es = SpatialESArray(wfArray)
np.save(fe % saveParas, ee)
np.save(fs % saveParas, es)
