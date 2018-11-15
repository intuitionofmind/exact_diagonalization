#!/usr/bin/python3

import os
import math
import cmath as cm
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from tjsquare import dim, LoadEigVal, LoadEigVec

fs = 12
inter = 5

def FullTimeEvo(t):
    U = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        lam = w[i]-w[0]
        U += np.exp(-1.j*lam*t)*np.outer(v[:, i], v[:, i])
    return U

def CutoffTimeEvo(t, cf):
    U = np.zeros((dim, dim), dtype=complex)
    for i in range(cf):
        lam = w[i]-w[0]
        U += np.exp(-1.j*lam*t)*np.outer(v[:, i], v[:, i])
    return U

def Normalize(v):
    z = np.absolute(np.vdot(v, v))
    return v/np.sqrt(z)


size = 52
dataDir = '/Users/wayne/Downloads/data/'
eigValsFile = 'eigenvalues1000_size52_J0.1_0.2_250_bcOO_sigma%s.dat'
eigVecsFile = 'eigenvectors1000_size52_J0.1_0.2_250_bcOO_sigma%s.dat'

parasARPACK = 'False'
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % parasARPACK))
wfArray = LoadEigVec(os.path.join(dataDir, eigVecsFile % parasARPACK))

# rawData = np.loadtxt('hamiltonian_sigmaTrue.txt')
matrFile = 'data/HamMatrix_size52_J%s_bcOO_sigma%s.txt'
J = 0.3
paras = (J, 'False')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
w, v = la.eigh(H)
print(w[:10])
print(ene[1][:10])

'''
vInit = v0[:, 1]
for i in range(numT):
    t = deltaT*(i+1)
    U = FullTimeEvo(t)
    V = CutoffTimeEvo(t, 200)
    matrixU = np.matrix(U)
    # print(i)
    # print(np.dot(matrixU, matrixU.getH()))
    vFinal = np.dot(V, vInit)
    # vFinal = Normalize(vFinal)
    # print(i, np.vdot(vInit, vFinal), np.vdot(vFinal, vFinal))
    print(i, np.vdot(np.dot(U, vInit), np.dot(V, vInit)))
'''
