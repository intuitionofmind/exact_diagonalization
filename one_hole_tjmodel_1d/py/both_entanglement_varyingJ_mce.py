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
from tjsquare_entanglement import BipartiteEE, MutualEE

size = 44
m = 100

dataDir = '/Users/wayne/Downloads/data/'

eigValsFile = 'eigenvalues100_size62_J0.1_0.1_100_bcOO_sigma%s.dat'
eigVecsFile = 'eigenvectors100_size62_J0.1_0.1_100_bcOO_sigma%s.dat'

# paras = (numEval, size, numSam, 'PP', 'False')
paras = 'False'
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile % paras))

ns = 100
size = 62
csvParas = (size, 'OO', paras)
fn = ['J', 'bee', 'mee']
writer = csv.DictWriter(open('both_ee_varyingJ_size%s_bc%s_sigma%s.csv' % csvParas, 'w'), fn)
writer.writeheader()

targetE = -5.0
deltaE = 0.02
J = 0.1
for i in range(ns):
    Jc = '{:.2f}'.format(J)
    be = 0.0
    me = 0.0
    n = 0
    for j in range(numSam):
        e = ene[i][j]
        if (targetE-deltaE)*J < e < (targetE+deltaE)*J:
            wf = wfArray[i][j]
            be += BipartiteEE(wf)
            me += MutualEE(wf)
            n += 1
    print(n)
    be = be/n
    me = me/n
    writer.writerow({'J': Jc, 'bee': be, 'mee': me})
    J = J+0.1
