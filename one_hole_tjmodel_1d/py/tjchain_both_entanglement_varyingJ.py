#!/usr/bin/python3

import os
import math
import random
import csv
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from tjchain import numSite, subDim, dim, holeBasis, spinBasis
from tjchain import FindState, LoadEigVal, LoadEigVec, Translation
from tjchain_entanglement import SpatialEE, MutualEE

m = 100

dataDir = '/Users/wayne/Downloads/data/'

eigValsFile = 'eigenvalues10_size1D10_J0.1_0.1_400_PBC_sigma%s.dat'
eigVecsFile = 'eigenvectors10_size1D10_J0.1_0.1_400_PBC_sigma%s.dat'

# paras = (numEval, size, numSam, 'PP', 'False')
paras = 'False'
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile % paras))

ns = 200
size = '1D10'
csvParas = (size, 'P', paras)
fn = ['J', 'bee', 'mee']
writer = csv.DictWriter(open('both_ee_varyingJ_size%s_bc%s_sigma%s.csv' % csvParas, 'w'), fn)
writer.writeheader()

J = 0.1
for i in range(ns):
    print(i)
    Jc = '{:.2f}'.format(J)
    wfList = [wfArray[i][0], wfArray[i][1]]
    fold = 2
    wfListNew = np.zeros((fold, dim), dtype=complex)
    # print(ene[i][0], ene[i][1])
    if np.absolute(ene[i][1]-ene[i][0]) < 1e-10:
        T = np.zeros((fold, fold), dtype=complex)
        for i in range(fold):
            for j in range(fold):
                T[i][j] = np.vdot(wfList[i], Translation(wfList[j]))
        U, X = la.schur(T)

        mom = []
        for i in range(fold):
            mom.append(np.angle(U[i][i]))

        for i in range(fold):
            for j in range(fold):
                wfListNew[i] += X[:, i][j]*wfList[j]

        # print(mom)
        if mom[0] > 0.0:
            wf = wfListNew[0]
        else:
            wf = wfListNew[1]
    else:
        wf = wfList[0]

    # be = BipartiteEE(wf)
    be = SpatialEE(wf)
    me = MutualEE(wf)
    print(i, be, me)
    writer.writerow({'J': Jc, 'bee': be, 'mee': me})
    J = J+0.1
