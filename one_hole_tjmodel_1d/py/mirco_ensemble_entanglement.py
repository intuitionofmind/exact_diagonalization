#!/usr/bin/python3

import os
import math
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from tjchain import dim, cut, numSite, numEval, numSam
from tjchain import LoadEigVal, LoadEigVec

dataDir = '/Users/wayne/Downloads/data/'
eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'
eigVecsFile = 'eigenvectors%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'

matrFile = 'data/HamMatrix_size1D10_J%s_bcP_sigma%s.txt'

level = 49

meeFile = 'mutualEE%s_len%s_J%s_%s_sigma%s.npy'
mesFile = 'mutualES%s_len%s_J%s_%s_sigma%s.npy'
seeFile = 'spatialEE%s_len%s_cut%s_J%s_%s_sigma%s.npy'
sesFile = 'spatialES%s_len%s_cut%s_J%s_%s_sigma%s.npy'

def MeasureSEE(paras, eneSpectrum, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    see = np.load(os.path.join(dataDir, seeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[i])/numSite
        if  (targetE-deltaE) < e < (targetE+deltaE): 
            temp = np.append(temp, see[i])
    # print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))

def MeasureMEE(paras, eneSpectrum, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    mee = np.load(os.path.join(dataDir, meeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[i])/numSite
        if (targetE-deltaE) < e < (targetE+deltaE):
            temp = np.append(temp, mee[i])
    # print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))


ParasF = (numEval, numSite, numSam, 'PBC', 'False')
eF = os.path.join(dataDir, eigValsFile) % ParasF
fF = os.path.join(dataDir, eigVecsFile) % ParasF
eneArrayF = LoadEigVal(eF)
wfArrayF = LoadEigVec(fF)

ParasT = (numEval, numSite, numSam, 'PBC', 'True')
eT = os.path.join(dataDir, eigValsFile) % ParasT
fT = os.path.join(dataDir, eigVecsFile) % ParasT
eneArrayT = LoadEigVal(eT)
wfArrayT = LoadEigVec(fT)

J = 0.3
deltaE = 0.01 # Resolution of the energy density

paras = (J, 'F')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
vInit = wfArrayF[4][level]
targetEFalse = np.vdot(vInit, np.dot(H, vInit))/numSite

mEEnpyParas = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
eneSpectrum = eneArrayF[0]

print('targetEFalse', targetEFalse)
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetEFalse, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetEFalse, deltaE))

paras = (J, 'T')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
vInit = wfArrayT[4][level]
targetETrue = np.vdot(vInit, np.dot(H, vInit))/numSite

mEEnpyParas = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
eneSpectrum = eneArrayT[0]

print('targetETrue', targetETrue)
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetETrue, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetETrue, deltaE))


J = 40.3 
deltaE = 0.4 # Resolution of the energy density

paras = (J, 'F')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
vInit = wfArrayF[0][level]
targetEFalse = np.vdot(vInit, np.dot(H, vInit))/numSite

mEEnpyParas = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
eneSpectrum = eneArrayF[4]

print('targetEFalse', targetEFalse)
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetEFalse, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetEFalse, deltaE))

paras = (J, 'T')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
vInit = wfArrayT[0][level]
targetETrue = np.vdot(vInit, np.dot(H, vInit))/numSite

mEEnpyParas = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
eneSpectrum = eneArrayT[4]

print('targetETrue', targetETrue)
print(sEEnpyParas)
print(MeasureMEE(mEEnpyParas, eneSpectrum, targetETrue, deltaE))
print(MeasureSEE(sEEnpyParas, eneSpectrum, targetETrue, deltaE))
