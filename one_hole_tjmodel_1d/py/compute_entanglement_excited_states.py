#!/usr/bin/python3

import os
import math
import numpy as np
from scipy import linalg as la

from auxi import numSite, numEval, cut, dim

from tjchain_entanglement import SpatialEEArray
from tjchain_entanglement import SpatialESArray

from tjchain_entanglement import MutualEEArray
from tjchain_entanglement import MutualESArray

fv = '/Users/wayne/Downloads/data/eigenvectors%s_size%s_J0.3_10.0_5_%s_sigma%s.dat'

def LoadEigVec(paras, l):
    rawData = np.fromfile(fv % paras, dtype=np.float)
    eigVec = np.zeros((paras[0], dim), dtype=complex)
    for i in range(paras[0]):
        for j in range(dim):
            eigVec[i][j] = complex(rawData[l*(paras[0]*dim*2)+i*dim*2+j*2], rawData[l*(paras[0]*dim*2)+i*dim*2+j*2+1])
    return eigVec

# Set the slice in the array.
l = 4
J = 40.3

# Compute spatial bibpartite entanglement.

fileSpatialEE = 'spatialEE%s_len%s_cut%s_J%s_%s_sigma%s'
fileSpatialES = 'spatialES%s_len%s_cut%s_J%s_%s_sigma%s'
fileMutualEE = 'mutualEE%s_len%s_J%s_%s_sigma%s'
fileMutualES = 'mutualES%s_len%s_J%s_%s_sigma%s'

loadParas = (numEval, numSite, 'PBC', 'False')
wfArray = LoadEigVec(loadParas, l)
spatialEE = SpatialEEArray(wfArray)
# spatialES = SpatialESArray(wfArray)
mutualEE = MutualEEArray(wfArray)
# mutualES = MutualESArray(wfArray)

saveSpatial = (numEval, numSite, cut, J, 'PBC', 'False')
saveMutual = (numEval, numSite, J, 'PBC', 'False')
np.save(fileSpatialEE % saveSpatial, spatialEE)
# np.save(fileSpatialES % saveSpatial, spatialES)
np.save(fileMutualEE % saveMutual, mutualEE)
# np.save(fileMutualES % saveMutual, mutualES)

# for sigma t-J model

loadParas = (numEval, numSite, 'PBC', 'True')
wfArray = LoadEigVec(loadParas, l)
spatialEE = SpatialEEArray(wfArray)
# spatialES = SpatialESArray(wfArray)
mutualEE = MutualEEArray(wfArray)
# mutualES = MutualESArray(wfArray)

saveSpatial = (numEval, numSite, cut, J, 'PBC', 'True')
saveMutual = (numEval, numSite, J, 'PBC', 'True')
np.save(fileSpatialEE % saveSpatial, spatialEE)
# np.save(fileSpatialES % saveSpatial, spatialES)
np.save(fileMutualEE % saveMutual, mutualEE)
# np.save(fileMutualES % saveMutual, mutualES)
