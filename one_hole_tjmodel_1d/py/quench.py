#!/usr/bin/python3

import os
import math
import cmath as cm
import numpy as np
from scipy import linalg as la
import csv

from auxi_1d import dim
from tjchain_entanglement import MutualEE, SpatialEE


def FullTimeEvo(w, v, t):
    U = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        lam = w[i]-w[0]
        U += np.exp(-1.j*lam*t)*np.outer(v[:, i], v[:, i])
    return U

def CutoffTimeEvo(w, v, t, cf):
    U = np.zeros((dim, dim), dtype=complex)
    for i in range(cf):
        lam = w[i]-w[0]
        U += np.exp(-1.j*lam*t)*np.outer(v[:, i], v[:, i])
    return U

def Normalize(v):
    z = np.absolute(np.vdot(v, v))
    return v/np.sqrt(z)

matrFile = 'data/HamMatrix_size1D10_J%s_bcP_sigma%s.txt'

JQ = 40.3
parasQ = (JQ, 'F')
rawDataQ = np.loadtxt(matrFile % parasQ)
HQ = np.reshape(rawDataQ, (dim, dim))
wQ, vQ = la.eigh(HQ)
print(wQ)

J = 0.3
paras = (J, 'F')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
w, v = la.eigh(H)
print(w)

deltaT = 0.1
numT = 1000
level = 100
vInit = vQ[:, level]

csvParas = (deltaT, level, 'P', 'F')
fn = ['num', 'mee', 'see']
writer = csv.DictWriter(open('timeEE_J40.3to0.3_step%s_state%s_bc%s_sigma%s.csv' % csvParas, 'w'), fn)
writer.writeheader()

for i in range(numT):
    t = deltaT*i
    U = FullTimeEvo(w, v, t)
    # V = CuitoffTimeEvo(w, v, t, 1000)
    # matrixU = np.matrix(U)
    # print(np.dot(matrixU, matrixU.getH()))
    vFinal = np.dot(U, vInit)
    mee = MutualEE(vFinal)
    see = SpatialEE(vFinal)
    print(i, mee, see)
    # vFinal= Normalize(vFinal)
    # print(i, np.vdot(vInit, np.dot(U, vInit)), np.vdot(vInit, np.dot(V, vInit)))
    writer.writerow({'num': i, 'mee': np.real(mee), 'see': np.real(see)})
