#!/usr/bin/python3

import os
import math
import cmath as cm
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

from auxi_1d import dim, subDim, numSite, holeBasis, spinBasis
fs = 12
inter = 5

def Sz(r, vec):
    w = np.zeros(dim, dtype=complex)
    for l in range(dim):
        h = holeBasis[l // subDim]
        s = spinBasis[l % subDim] #i = h*dim_sub+s
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'2'+b[h:] # To restore the t-J configuration. 
        if config[r] == '0':
            w[l] = -0.5*vec[l]
        elif config[r] == '1':
            w[l] = 0.5*vec[l]
    return w

def FullTimeEvo(t, w, v):
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

def SzCommutatorSquare(r0, r1, t, vec):
   # w(t)vvw(t)
    p0 = np.vdot(vec, np.dot(UInver, Sz(r0, np.dot(U, Sz(r1, Sz(r1, np.dot(UInver, Sz(r0, np.dot(U, vec)))))))))
    # vw(t)w(t)v
    p1 = np.vdot(vec, Sz(r1, np.dot(UInver, Sz(r0, Sz(r0, np.dot(U, Sz(r1, vec)))))))
    # vw(t)vw(t)
    p2 = np.vdot(vec, Sz(r1, np.dot(UInver, Sz(r0, np.dot(U, Sz(r1, np.dot(UInver, Sz(r0, np.dot(U, vec)))))))))
    # w(t)vw(t)v
    p3 = np.vdot(vec, np.dot(UInver, Sz(r0, np.dot(U, Sz(r1, np.dot(UInver, Sz(r0, np.dot(U, Sz(r1, vec)))))))))
    return p0+p1-p2-p3


matrFile = 'data/HamMatrix_size1D10_J%s_bcP_sigma%s.txt'
J = 0.3
paras = (J, 'F')
rawData = np.loadtxt(matrFile % paras)
H = np.reshape(rawData, (dim, dim))
eigVal, eigVec = la.eigh(H)
print(eigVal[:10])

beta = 20.0 # Inversed temperature. 
cut = int(dim/2)
num = 100
timeStep = 0.1

'''
x0 = 1
y0 = 0
r0 = y0*numSiteX+x0
x1 = 2
y1 = 1
r1 = y1*numSiteX+x1
'''

r0 = 2
r1 = 7

com = np.zeros(num, dtype=float)
for i in range(num):
    t = i*timeStep
    U = FullTimeEvo(t, eigVal, eigVec)
    UInver = FullTimeEvo(-1.*t, eigVal, eigVec)
    tr = 0.0
    pf = 0.0
    for j in range(cut):
        deltaE = eigVal[j]-eigVal[0]
        factor = np.exp(-1.0*beta*deltaE)
        v = eigVec[:, j]
        butterfly = SzCommutatorSquare(r0, r1, t, v)
        tr += factor*butterfly
        pf += factor
    tr = tr/pf
    com[i] = np.real(tr)
    print(t, tr)

f = 'otoc_sz_J%s_beta%s_sigma%s_num%s_step%s'
saveParas = (J, beta, sigma, num, timeStep)
np.save(f % saveParas, com)
