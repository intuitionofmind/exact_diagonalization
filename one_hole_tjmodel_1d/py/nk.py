#!/usr/bin/python3

import os
import math
import cmath as cm
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from auxi import *

# Hopping() computes the c_{i\sigma}^{\dagger}c_{j\sigma}|v\rangle.

def Hopping(v, sigma, i, j):
    w = np.zeros(dim, dtype=complex)
    for k in range(dim):
        h = holeBasis[k // subDim]
        s = spinBasis[k % subDim]
        b = bin(s)[2:]
        b = '%0*.0f' % (numSite, int(b)) 
        b = b[::-1]
        config = b[:h]+sigma+b[h:] # temporary half filled config 
        if i == j:
            if i != h and sigma == config[j]:
                w[k] += v[k]
        else:
            if i == h and sigma == config[j]:
                hh = j
                bb = config[:hh]+config[hh+1:]
                bb = bb[::-1]
                ss = int(bb, 2)
                l = FindState(ss, 0, subDim-1)
                w[hh*subDim+l] += np.power(-1, i+j+1)*v[k]
    return w

def HoppingMatrix(v, sigma):
    M = np.zeros((numSite, numSite), dtype=complex)
    for i in range(numSite):
        for j in range(numSite):
            M[i][j] = np.vdot(v, Hopping(v, sigma, i, j))
    # print('hopping matrix', sigma, np.trace(M))
    return M

def ProjectedHole(v, hh):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        h = holeBasis[i // subDim]
        if h != hh:
            continue
        else:
            w[i] = v[i]
    return w

def HoleDistribution(paras, i):
    T = np.zeros((2, 2), dtype=complex)
    m = np.zeros(numSite, dtype=complex)
    if np.absolute(ene[i][0]-ene[i][1]) < 1e-10:
        wf0t = Translation(wf0[i])
        wf1t = Translation(wf1[i])
        Tem[0][0] = np.vdot(wf0[i], wf0t) # vdot for complex conjugate...
        Tem[0][1] = np.vdot(wf0[i], wf1t)
        Tem[1][0] = np.vdot(wf1[i], wf0t)
        Tem[1][1] = np.vdot(wf1[i], wf1t)
        w, v = la.eig(T)
        print(i, np.angle(w, deg=True))
        if w[0].imag > 0:
            choice = 1
        else:
            choice = 0
        wf = v[:, choice][0]*wf0[i]+v[:, choice][1]*wf1[i]  # construct eigenstates with specific momentum
    else:
        wf = wf0[i]
    for l in range(numSite):
        m[l] = np.vdot(wf, ProjectedHole(wf, l))
    return np.real(m)

def NkSample(i, sigma):
    T = np.zeros((2, 2), dtype=complex)
    m = np.zeros(numSite, dtype=complex)
    if np.absolute(ene[i][0]-ene[i][1]) < 1e-10:
        wf0t = Translation(wf0[i])
        wf1t = Translation(wf1[i])
        T[0][0] = np.vdot(wf0[i], wf0t) # vdot for complex conjugate...
        T[0][1] = np.vdot(wf0[i], wf1t)
        T[1][0] = np.vdot(wf1[i], wf0t)
        T[1][1] = np.vdot(wf1[i], wf1t)
        w, v = la.eig(T)
        print(i, np.angle(w, deg=True))
        if w[0].imag > 0: # choose the eige
            choice = 0
        else:
            choice = 1
        wf = v[:, choice][0]*wf0[i]+v[:, choice][1]*wf1[i]  # construct eigenstates with specific momentum
    else:
        wf = wf0[i]

    M = HoppingMatrix(wf, sigma)
    for l in range(numSite):
        for j in range(numSite):
            for k in range(numSite):
                m[l] += np.exp(1j*(2*np.pi*l/numSite)*(j-k))*M[j][k]/numSite
    print(sigma, np.sum(m))
    return np.real(m)

# --------------Test begin--------------------

# # PrintHam(1.1, False)
# slice = 10
# paras = (numEval, numSite, numSam, 'PBC', 'false', flux)
# print(paras)
# ene = LoadEigVal(paras)
# wf = LoadEigVec(paras, 0)
# v = wf[slice]
# print('Check energy: ', ene[slice][0], np.vdot(v, Hamiltonian(v, (slice+1)*0.1, False)))
# # print(HoleDistribution(paras, slice))

#---------------Test end-----------------------

fig = plt.figure(figsize=(8, 4))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
x = np.arange(numSite)*2*np.pi/numSite
u = np.ones(numSite, dtype=float)
mark = [1.0]*numSite

slice1 = 350
slice2 = 10
paras = (numEval, numSite, numSam, 'PBC', 'false', flux)
print(paras)
ene = LoadEigVal(paras)
wf0 = LoadEigVec(paras, 0)
wf1 = LoadEigVec(paras, 1)
m1 = NkSample(slice1, '1')+NkSample(slice1, '0')
m2 = NkSample(slice2, '1')+NkSample(slice2, '0')

paras = (numEval, numSite, numSam, 'PBC', 'true', flux)
print(paras)
ene = LoadEigVal(paras)
wf0 = LoadEigVec(paras, 0)
wf1 = LoadEigVec(paras, 1)
mm1 = NkSample(slice1, '1')+NkSample(slice1, '0')
mm2 = NkSample(slice2, '1')+NkSample(slice2, '0')

ax1 = plt.subplot(121)
plt.plot(x, m1, '-o', color='red', label='$t$-$J$')
plt.plot(x, mm1, '-->', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.axhline(y=1.0, xmin=0.0, xmax=1.0, linestyle='-.', color='black')
plt.legend(loc='best', frameon=False, fontsize=fs)
# plt.setp(ax1.get_xticklabels(), visible=False)
plt.xticks([0, np.pi, 2*np.pi],
        ['$0$', '$\pi$', '$2\pi$'], fontsize=fs)
plt.xlabel('$k$', fontsize=fs)
plt.ylabel('$n_{k}$', fontsize=20)
# start, end = ax1.get_ylim()
start = 0.3
end = 1.3
ax1.set_ylim([start, end])
ax1.set_yticks(np.arange(start, end, 0.2))
plt.setp(ax1.get_yticklabels(), fontsize=fs)
ax1.set_title('(a) $J/t=35.0$', fontsize=fs, x=0.25, y=0.85)

# ax2 = plt.subplot(222, sharey=ax1)
# m = m1-mm1
# print(m)
# plt.plot(x, m, '-o', color='red', label='$n_{k}^{t-J}-n_{k}^{\sigma\cdot{t}-J}, J/t=1.0$')
# plt.plot(x, mark, '--', color='blue')
# plt.legend(loc='best', frameon=False, fontsize=12)
# plt.setp(ax2.get_xticklabels(), visible=False)
# plt.setp(ax2.get_yticklabels(), visible=False)
# plt.ylabel('$n_{k}^{t-J}-n_{k}^{\sigma\cdot{t}-J}$', fontsize=14)
# ax2.set_title('(b)', fontsize=14)

ax3 = plt.subplot(122, sharey=ax1)
plt.plot(x, m2, '-o', color='red', label='$t$-$J$')
plt.plot(x, mm2, '-->', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.axhline(y=1.0, xmin=0.0, xmax=1.0, linestyle='-.', color='black')
plt.legend(loc='best', frameon=False, fontsize=fs)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.xticks([0, np.pi, 2*np.pi],
        ['$0$', '$\pi$', '$2\pi$'], fontsize=fs)
plt.xlabel('$k$', fontsize=fs)
# plt.ylabel('$n_{k}$', fontsize=fs)
start = 0.3
end = 1.3
ax3.set_ylim([start, end])
ax3.set_yticks(np.arange(start, end, 0.2))
plt.setp(ax3.get_yticklabels(), fontsize=fs)
ax3.set_title('(b) $J/t=1.0$',fontsize=fs, x=0.25, y=0.85)

# ax4 = plt.subplot(224, sharex=ax2, sharey=ax3)
# m = m2-mm2
# print(m)
# plt.plot(x, m, '-o', color='red', label='$n_{k}^{t-J}-n_{k}^{\sigma\cdot{t}-J}, J/t=35.0$')
# plt.plot(x, mark, '--', color='blue')
# plt.legend(loc='best', frameon=False, fontsize=12)
# plt.setp(ax4.get_yticklabels(), visible=False)
# plt.xticks([0, np.pi, 2*np.pi],
        # ['$0$', '$\pi$', '$2\pi$'], fontsize=14)
# plt.xlabel('$k$', fontsize=14)
# plt.ylabel('$n_{k}^{t-J}-n_{k}^{\sigma\cdot{t}-J}$', fontsize=14)
# start = -0.1
# end = 1.1
# ax3.set_ylim([start, end])
# ax4.set_title('(d)')

image = 'nk_len_%s_slice_%s_flux_%s.pdf'
imageParas = (numSite, slice1, flux)
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
plt.show()
