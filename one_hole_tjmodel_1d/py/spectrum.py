#!/usr/bin/python3

import os
import math
import numpy as np

from tjchain import numSite

dataDir = '/Users/wayne/Downloads/data/'
eigvalsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'

def LoadEigvalues(f):
    rawData = np.fromfile(f, dtype=np.float)
    return np.reshape(rawData, (numSam, numEval))

numEval = 1200
numSam = 5
paras = (numEval, numSite, numSam, 'PBC', 'False')
# paras = (numEval, numSite, numSam, 'PBC', 'False')
f = os.path.join(dataDir, eigvalsFile) % paras
sp = LoadEigvalues(f)

l = 99
print(sp[0][l])


