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

dataDir = '/Users/wayne/Downloads/data/'

eigValsFile = 'eigenvalues1000_size62_J0.1_0.1_10_bcOO_sigma%s.dat'
eigVecsFile = 'eigenvectors1000_size62_J0.1_0.1_10_bcOO_sigma%s.dat'
# eigValsFile = 'eigenvalues1000_size44_J0.3_10.0_5_PP_sigma%s.dat'
# eigVecsFile = 'eigenvectors1000_size44_J0.3_10.0_5_PP_sigma%s.dat'



# paras = (numEval, size, numSam, 'PP', 'False')
paras = 'False'
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile % paras))

ns = 1000
csvParas = (size, 'OO', paras)
fn = ['ene_density', 'bee', 'mee']
writer = csv.DictWriter(open('both_ee_excited_size%s_bc%s_sigma%s.csv' % csvParas, 'w'), fn)
writer.writeheader()

l = 9
for i in range(ns):
    e = ene[l][i]/18.0
    wf = wfArray[l][i]
    be = BipartiteEE(wf)
    me = MutualEE(wf)
    writer.writerow({'ene_density': e, 'bee': be, 'mee': me})
