#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from auxi import *

def RedDenMatrixHole(wf):
    Rho = np.zeros((numSite, numSite), dtype=complex)
    for i in range(0, numSite):
        for j in range(0, numSite):
            for k in range(0, subDim):
                Rho[i][j] += wf[i*subDim+k]*np.conjugate(wf[j*subDim+k])
    return Rho

def MutualES(wf):
    Rho = RedDenMatrixHole(wf)
    es = np.linalg.eigvalsh(Rho)
    return es

def MutualEE(wf):
    es = MutualES(wf)
    ee = -1.0*np.dot(es, np.log(np.absolute(es)))  # entanglement entropy
    return ee

def MutualEEArray(wfArray):
    eeArray = np.zeros(numEval)
    for i in range(numEval):
        eeArray[i] = MutualEE(wfArray[i])
        print(i, eeArray[i])
    return eeArray

def MutualESArray(wfArray):
    esArray = np.zeros((numEval, numSite))
    for i in range(numEval):
        esArray[i] = MutualES(wfArray[i])
    return esArray 

l = 9
fe = 'mutual-ee_len-%s_J-%s_%s_sigma-%s'
fs = 'mutual-es_len-%s_J-%s_%s_sigma-%s'

loadParas = (numEval, numSite, numSam, 'PBC', 'false')
saveParas = (numSite, (l+1)*step, 'PBC', 'false')
wfArray = LoadEigVec2(loadParas, l)
ee = MutualEEArray(wfArray)
es = MutualESArray(wfArray)
np.save(fe % saveParas, ee)
np.save(fs % saveParas, es)

loadParas = (numEval, numSite, numSam, 'PBC', 'true')
saveParas = (numSite, (l+1)*step, 'PBC', 'true')
wfArray = LoadEigVec2(loadParas, l)
ee = MutualEEArray(wfArray)
es = MutualESArray(wfArray)
np.save(fe % saveParas, ee)
np.save(fs % saveParas, es)
