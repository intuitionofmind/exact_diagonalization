#!/usr/bin/python3

import os
import math
import numpy as np
import logging
from scipy import linalg as la

from tjchain import numSite, numEval, dim, subDim, spinBasis, holeBasis


# --------------------- bipartite entanglement -- start --------------------
# This section consists of functions to compute the bipartite entanglement for $t$-$J$ model in one dimension.

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
        print(config, configA, configB)
        if configA not in basisA:
            basisA.append(configA)
        if configB not in basisB:
            basisB.append(configB)
        indexDict.append((basisA.index(configA), basisB.index(configB)))

dimA = len(basisA)
dimB = len(basisB)
print(dimA, dimB, dim)

def Reshape(wf):
    WF = np.zeros((dimA, dimB), dtype=complex)
    for i in range(dim):
        j = indexDict[i][0]
        k = indexDict[i][1]
        WF[j][k] = wf[i]
    return WF

def RedDenMatrixA(WF):
    WF = np.matrix(WF)
    return WF*(WF.getH())

def SpatialES(wf):
    WF = Reshape(wf)
    # r = np.outer(np.conjugate(wf), wf)
    # print(np.linalg.matrix_rank(r), np.linalg.matrix_rank(WF))
    Rho = RedDenMatrixA(WF)
    es = np.linalg.eigvalsh(Rho)
    l = len(es)
    # Set a small but finite trunction for the sake of the $log$ computation.
    for i in range(l):
        if 0.0 == es[i]:
            es[i] = 1e-99
    return es

def SpatialEH(wf):
    es = SpatialES(wf)
    return (-1.0*np.log(es))

def SpatialEE(wf):
    es = SpatialES(wf)
    ee = -1.0*np.dot(es, np.log(np.absolute(es)))  # entanglement entropy
    # print(ee)
    return ee

def SpatialEEArray(wfArray):
    eeArray = np.zeros(numEval)
    for i in range(numEval):
        eeArray[i] = SpatialEE(wfArray[i])
        print('sEE', i, eeArray[i])
    return eeArray

def SpatialESArray(wfArray):
    esArray = np.zeros((numEval, dimA))
    for i in range(numEval):
        esArray[i] = SpatialES(wfArray[i])
    return esArray 

# --------------------- bipartite entanglement -- end --------------------


# --------------------- mutual entanglement -- start --------------------

def ReshapeWaveFunctionHole(wf):
    WF = np.zeros((numSite, subDim), dtype=complex)
    for i in range(numSite):
        for j in range(subDim):
            k = i*subDim+j
            WF[i][j] = wf[k]
    return WF

def RedDenMatrixHole(wf):
    WF = ReshapeWaveFunctionHole(wf)
    WF = np.matrix(WF)
    return WF*(WF.getH())

def MutualES(wf):
    Rho = RedDenMatrixHole(wf)
    es = np.linalg.eigvalsh(Rho)
    l = len(es)
    # Set a small but finite trunction for the sake of the $log$ computation.
    for i in range(l):
        if 0.0 == es[i]:
            es[i] = 1e-99
    return es

def MutualEH(wf):
    es = MutualES(wf)
    return (-1.0*np.log(es))

def MutualEE(wf):
    es = MutualES(wf)
    ee = -1.0*np.dot(es, np.log(np.absolute(es)))  # entanglement entropy
    return ee

def MutualEEArray(wfArray):
    eeArray = np.zeros(numEval)
    for i in range(numEval):
        eeArray[i] = MutualEE(wfArray[i])
        # print('mEE', i, eeArray[i])
    return eeArray

def MutualESArray(wfArray):
    esArray = np.zeros((numEval, numSite))
    for i in range(numEval):
        esArray[i] = MutualES(wfArray[i])
    return esArray 

# --------------------- mutual entanglement -- end --------------------
