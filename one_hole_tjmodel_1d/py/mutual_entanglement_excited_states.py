#!/usr/bin/python3

import os
import math
import numpy as np
from scipy import linalg as la

from auxi import *

from tjchain_entanglement import 


l = 0
J = 1.0
fe = 'mutual_EE%s_len%s_J%s_%s_sigma%s'
fs = 'mutual_ES%s_len%s_J%s_%s_sigma%s'

loadParas = (numEval, numSite, numSam, 'PBC', 'false')
saveParas = (numEval, numSite, J, 'PBC', 'False')
wfArray = LoadEigVec2(loadParas, l)
ee = MutualEEArray(wfArray)
es = MutualESArray(wfArray)
np.save(fe % saveParas, ee)
np.save(fs % saveParas, es)

loadParas = (numEval, numSite, numSam, 'PBC', 'true')
saveParas = (numEval, numSite, J, 'PBC', 'True')
wfArray = LoadEigVec2(loadParas, l)
ee = SpatialEEArray(wfArray)
es = SpatialESArray(wfArray)
np.save(fe % saveParas, ee)
np.save(fs % saveParas, es)
