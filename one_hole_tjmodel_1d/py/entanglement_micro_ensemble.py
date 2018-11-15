#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid.inset_locator import inset_axes
from auxi import LoadEigVal
from auxi import Hamiltonian
from auxi import cut

fs = 16 # font size 

dataDir = '/Users/wayne/Downloads/data/'

numSite = 10
numEval = 1200
level = 220
numSam = 5

J = 0.3
l = 0
eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'

meeFile = 'mutualEE%s_len%s_J%s_%s_sigma%s.npy'
mesFile = 'mutualES%s_len%s_J%s_%s_sigma%s.npy'
seeFile = 'spatialEE%s_len%s_cut%s_J%s_%s_sigma%s.npy'
sesFile = 'spatialES%s_len%s_cut%s_J%s_%s_sigma%s.npy'

def MeasureSEE(paras, eneSpectrum, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    see = np.load(os.path.join(dataDir, seeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[l][i])/numSite
        if  (targetE-deltaE) < e < (targetE+deltaE): 
            temp = np.append(temp, see[i])
    # print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))

def MeasureMEE(paras, eneSpectrum, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    mee = np.load(os.path.join(dataDir, meeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[l][i])/numSite
        if (targetE-deltaE) < e < (targetE+deltaE):
            temp = np.append(temp, mee[i])
    # print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))

eneParas = (numEval, numSite, numSam, 'PBC', 'False')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), eneParas)
# print(eneSpectrum)
mEEnpyParas = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
targetEFalse = (-1.7816)/numSite
deltaE = 0.01 # Resolution of the energy density
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetEFalse, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetEFalse, deltaE))

eneParas = (numEval, numSite, numSam, 'PBC', 'True')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), eneParas)
# print(eneSpectrum)
mEEnpyParas = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
targetETrue = (-1.7008)/numSite
deltaE = 0.01 # Resolution of the energy density
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetETrue, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetETrue, deltaE))

J = 40.3 
l = 4

eneParas = (numEval, numSite, numSam, 'PBC', 'False')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), eneParas)
# print(eneSpectrum[l])
mEEnpyParas = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
targetEFalse = (-79.05)/numSite
deltaE = 1.0 # Resolution of the energy density
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetEFalse, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetEFalse, deltaE))

eneParas = (numEval, numSite, numSam, 'PBC', 'True')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), eneParas)
mEEnpyParas = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
targetETrue = (-36.01)/numSite
deltaE = 1.0 # Resolution of the energy density
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetETrue, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetETrue, deltaE))

