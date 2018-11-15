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

size = 62
m = 100

dataDir = '/Users/wayne/Downloads/data/'

eigValsFile = 'eigenvalues100_size62_J0.1_0.1_100_bcOO_sigma%s.dat'
eigVecsFile = 'eigenvectors100_size62_J0.1_0.1_100_bcOO_sigma%s.dat'

# paras = (numEval, size, numSam, 'PP', 'False')
paras = 'True'
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile % paras))

ns = 100
size = 62
csvParas = (size, 'OO', paras)
fn = ['J', 'bee', 'mee']
writer = csv.DictWriter(open('both_ee_varyingJ_size%s_bc%s_sigma%s.csv' % csvParas, 'w'), fn)
writer.writeheader()

J = 0.1
for i in range(ns):
    Jc = '{:.2f}'.format(J)
    wf = wfArray[i][0]
    be = BipartiteEE(wf)
    me = MutualEE(wf)
    writer.writerow({'J': Jc, 'bee': be, 'mee': me})
    J = J+0.1
