#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la

from auxi_1d import numSite

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
