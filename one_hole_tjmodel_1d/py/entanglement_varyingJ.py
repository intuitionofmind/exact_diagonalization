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
from tjchain_entanglement import SpatialEE, MutualEE, SpatialEH, MutualEH


dataDir = '/Users/wayne/Downloads/data/'

eigValsFile = 'eigenvalues10_size1D10_J0.1_0.1_400_PBC_sigma%s.dat'
eigVecsFile = 'eigenvectors10_size1D10_J0.1_0.1_400_PBC_sigma%s.dat'

# paras = (numEval, size, numSam, 'PP', 'False')
paras = 'True'
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile % paras))

fold = 2
size = 10 
csvParas = (size, paras)
# fn = ['J', 'bee', 'mee']
fn = ['J', 'beh', 'meh']
writer = csv.DictWriter(open('eh_varyingJ_size%s_PBC_sigma%s.csv' % csvParas, 'w'), fn)
writer.writeheader()

J = 0.1
for i in range(400):
    Jc = '{:.2f}'.format(J)
    if (ene[i][1]-ene[i][0]) < 1e-10:
        wfArraySub = wfArray[i][:fold]
        T = np.zeros((fold, fold), dtype=complex)
        for i in range(fold):
            for j in range(fold):
                T[i][j] = np.vdot(wfArraySub[i], Translation(wfArraySub[j]))
        U, X = la.schur(T) # T = XUX^{\dagger} 
        # print(Jc, np.angle(U[0][0], deg=True), np.angle(U[1][1], deg=True))
        wfArrayTrans = np.zeros((fold, dim), dtype=complex)
        for i in range(fold):
            for j in range(fold):
                wfArrayTrans[i] += X[:, i][j]*wfArraySub[j]
        if np.angle(U[0][0] > 0):
            wf = wfArrayTrans[0]
        else:
            wf = wfArrayTrans[1]

    else:
        wf = wfArray[i][0]
        phase = np.vdot(wf, Translation(wf))
        # print(Jc, np.angle(phase))
    be = SpatialEE(wf)
    me = MutualEE(wf)
    bh = SpatialEH(wf)
    mh = MutualEH(wf)
    # mh = MutualEH(wf)
    # print(Jc, mh)
    # print(Jc, be, me)
    # me = MutualEE(wf)
    # writer.writerow({'J': Jc, 'bee': be, 'mee': me})
    writer.writerow({'J': Jc, 'beh': bh[-1], 'meh': mh[-1]})
    J = J+0.1
