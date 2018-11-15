#!/usr/bin/python3

import os
import math
import random
import csv
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from tjsquare import numSite, subDim, dim, holeBasis, spinBasis
from tjsquare import FindState, LoadEigVal, LoadEigVec
from tjsquare_entanglement import BipartiteEE, MutualEE, ReshapeWaveFunction

def SwapColumns(myArray, col0, col1):
    temp = myArray[:, col0]
    myArray[:, col0] = myArray[:, col1]
    myArray[:, col1] = temp

size = 62

dataDir = '/Users/wayne/Downloads/data/'

eigValsFile = 'eigenvalues1000_size62_J0.1_0.1_10_bcOO_sigma%s.dat'
eigVecsFile = 'eigenvectors1000_size62_J0.1_0.1_10_bcOO_sigma%s.dat'

# paras = (numEval, size, numSam, 'PP', 'False')
paras = 'False'
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile % paras))

# e = ene[l][i]/18.0
wf = wfArray[0][0]
M = ReshapeWaveFunction(wf)
Rho0 = M*M.H
# print(Rh0)
es0 = np.linalg.eigvalsh(Rho0)
ee0 = -1.0*np.dot(es0, np.log(np.absolute(es0)))

# Rho1 = np.copy(Rho0)
N = np.copy(M)
N[:, [0, 8]] = N[:, [8, 0]]
N = np.matrix(N)
Rho1 = N*N.H
# print(Rho1)
# Rho1[:, [0, 1]] = Rho1[:, [1, 0]]
es1 = np.linalg.eigvalsh(Rho1)
ee1 = -1.0*np.dot(es1, np.log(np.absolute(es1)))

print(es0-es1)
print(ee0, ee1)
