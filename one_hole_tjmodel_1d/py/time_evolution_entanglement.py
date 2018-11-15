#!/usr/bin/python3

import os
import math
import cmath as cm
import numpy as np
import scipy as scp
from scipy import linalg as la

from auxi import Hamiltonian, LoadEigVal
from auxi import numSite, numEval
from tjchain_entanglement import MutualEE
from tjchain_entanglement import SpatialEE

import logging
logger = logging.getLogger('time_evolution')
hdlr = logging.FileHandler('time_evolution.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
logger.setLevel(logging.INFO)

dataDir = '/Users/wayne/Downloads/data/'

from auxi import numSite, dim, numSam
deltaTime = 0.001
num = 400
inter = 4
# inter = 2
loop = int(num/inter)

level = 240
l = 0
J = 40.3

eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'
fm = 'time_evolution_mEE%s_level%s_len%s_J%s_%s_sigma%s'
fs = 'time_evolution_sEE%s_level%s_len%s_J%s_%s_sigma%s'
statesFile = 'time_states_level%s_J%s_sigma%s.dat'

mee = np.zeros(loop, dtype=float)
see = np.zeros(loop, dtype=float)


sigma = False
paras = (numEval, numSite, numSam, 'PBC', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
print('ground state energy', ene[l][0])
saveParas = (num, level, numSite, J, 'PBC', 'False')
raw = np.fromfile(os.path.join(dataDir, statesFile) % (level, J, sigma))
wf = np.zeros(dim, dtype=complex)

for ii in range(int(num/inter)):
    i = ii*inter
    for j in range(dim):
        wf[j] = complex(raw[i*(dim*2)+j*2], raw[i*(dim*2)+j*2+1])
    mee[ii] = MutualEE(wf)
    see[ii] = SpatialEE(wf)
    print(sigma, i, wf[0], mee[ii], see[ii])

np.save(fm % saveParas, mee)
np.save(fs % saveParas, see)

# Turn to \sigma t-J model.

sigma = True 
paras = (numEval, numSite, numSam, 'PBC', 'True')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
print('ground state energy', ene[l][0])
saveParas = (num, level, numSite, J, 'PBC', 'True')
raw = np.fromfile(os.path.join(dataDir, statesFile) % (level, J, sigma))
wf = np.zeros(dim, dtype=complex)
for ii in range(loop):
    i = ii*inter
    for j in range(dim):
        wf[j] = complex(raw[i*(dim*2)+j*2], raw[i*(dim*2)+j*2+1])
    mee[ii] = MutualEE(wf)
    see[ii] = SpatialEE(wf)
    print(sigma, i, wf[0], mee[ii], see[ii])

np.save(fm % saveParas, mee)
np.save(fs % saveParas, see)
