#!/usr/bin/python3

import os
import math
import numpy as np
from scipy import linalg as la

from tjsquare import numSiteX, numSiteY, numSite, subDim, dim, holeBasis, spinBasis
from tjsquare import FindState, LoadEigVal, LoadEigVec

cut = 3

# Bipartite entanglement.

# Count the dimension of the Hilbert space for subsystems.
basisA = []
basisB = []
indexDict = [] # Construct the dictonary which maps indices between two representations. 

for l in range(dim):
    h = holeBasis[l // subDim] # h = hy*numSiteX+hx 
    hx = h % numSiteX
    hy = h // numSiteX
    s = spinBasis[l % subDim]
    b = bin(s)[2:]
    b = b.rjust(numSite-1, '0')
    b = b[::-1]
    config = b[:h]+'2'+b[h:]
    configA = ''
    configB = ''
    for j in range(numSiteY):
        # print(j, config, config[j*numSiteX:j*numSiteX+cut], config[j*numSiteX+cut:(j+1)*numSiteX])
        configA = configA+config[j*numSiteX:j*numSiteX+cut]
        configB = configB+config[j*numSiteX+cut:(j+1)*numSiteX]
    if configA not in basisA:
        basisA.append(configA)
    if configB not in basisB:
        basisB.append(configB)
    indexDict.append((basisA.index(configA), basisB.index(configB)))

dimA = len(basisA)
dimB = len(basisB)

print(dimA, dimB)

def ReshapeWaveFunction(wf):
    wfMatrix = np.zeros((dimA, dimB), dtype=complex)
    for l in range(dim):
        i = indexDict[l][0]
        j = indexDict[l][1]
        wfMatrix[i][j] = wf[l]
    return np.matrix(wfMatrix)

def BibpartiteES(wf):
    M = ReshapeWaveFunction(wf)
    Rho = M*M.H
    es = np.linalg.eigvalsh(Rho)
    for i in range(len(es)):
        if 0.0 == es[i]:
            es[i] = 1e-99
    return es

def BipartiteEE(wf):
    es = BibpartiteES(wf)
    ee = -1.0*np.dot(es, np.log(np.absolute(es)))
    return ee

# Mutual entanglement.

def MutualES(wf):
    Rho = np.zeros((numSite, numSite), dtype=complex)
    M = np.matrix(np.reshape(wf, (numSite, subDim))) # Reshape the wavefunction to the direct product matrix. 
    np.set_printoptions(precision=5, suppress=True)
    Rho = M*M.H # Matrix multiplication denotes the tracing out of the spin degrees of freedom. 
    print(np.real(Rho))
    es = np.linalg.eigvalsh(Rho)
    for i in range(len(es)):
        if 0.0 == es[i]:
            es[i] = 1e-99
    return es

def MutualEE(wf):
    es = MutualES(wf)
    ee = -1.0*np.dot(es, np.log(np.absolute(es)))  # entanglement entropy
    return ee

