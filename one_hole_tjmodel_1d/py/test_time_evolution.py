#!/usr/bin/python3

import os
import math
import cmath as cm
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from auxi_1d import dim
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

# rawData = np.loadtxt('hamiltonian_sigmaTrue.txt')
matrFile = 'data/HamMatrix_size1D10_J%s_bcP_sigma%s.txt'
J = 0.3
paras = (J, 'F')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
w, v = la.eigh(H)
print(w[:10])

'''
J = 10.3
paras = (J, 'OO', 'False')
# rawData = np.loadtxt('matrix_sigmaTrue.txt')
rawData = np.loadtxt(matrFile % paras)
# print(rawData)

H0 = np.zeros((dim, dim), dtype=complex)
for i in range(dim):
    for j in range(dim):
        H0[i][j] = complex(rawData[i*dim+j], 0.0)

w0, v0 = la.eigh(H0)
print(w0[:10])

deltaT = 0.01
numT = 1000

vInit = v0[:, 1]
for i in range(numT):
    t = deltaT*(i+1)
    U = FullTimeEvo(t)
    V = CutoffTimeEvo(t, 100)
    matrixU = np.matrix(U)
    # print(i)
    # print(np.dot(matrixU, matrixU.getH()))
    vFinal = np.dot(V, vInit)
    # vFinal = Normalize(vFinal)
    # print(i, np.vdot(vInit, vFinal), np.vdot(vFinal, vFinal))
    print(i, np.vdot(np.dot(U, vInit), np.dot(V, vInit)))
'''
