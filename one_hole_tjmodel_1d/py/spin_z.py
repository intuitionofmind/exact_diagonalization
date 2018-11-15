#!/usr/bin/python3

import os
import math
import cmath as cm
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from auxi import *

# Total Sz operator \sum_{i=0}^{L-1}S_{i}^{z}.

def TotalSz(v):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        s = spinBasis[i % subDim] #i = h*dim_sub+s
        h = holeBasis[i // subDim]
        b = bin(s)[2:]
        b = '%0*.0f' % (numSite, int(b))  #the highest bit is filled with a '0' just for convenience while does not have real meaning
        b = b[::-1]
        for j in range(numSite):
            if j < h:
                w[i] += (float(b[j])-0.5)*v[i]
            elif j > h:
                w[i] += (float(b[j-1])-0.5)*v[i]
            else:
                continue
    return w

# The z part of the square the total spin operator: \sum_{i=0}^{L-1}\sum_{j=0}^{L-1}S_{i}^{z}S_{j}^{z}.

def TotalSzSquare(v):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        s = spinBasis[i % subDim] #i = h*dim_sub+s
        h = holeBasis[i // subDim]
        b = bin(s)[2:]
        b = '%0*.0f' % (numSite, int(b))  #the highest bit is filled with a '0' just for convenience while does not have real meaning
        b = b[::-1]
        for j in range(numSite):
            for k in range(numSite):
                if j != h and k != h:
                    if j < h and k < h:
                        w[i] += (float(b[j])-0.5)*(float(b[k])-0.5)*v[i]
                    elif j < h and k > h:
                        w[i] += (float(b[j])-0.5)*(float(b[k-1])-0.5)*v[i]
                    elif j > h and k < h:
                        w[i] += (float(b[j-1])-0.5)*(float(b[k])-0.5)*v[i]
                    elif j > h and k > h:
                        w[i] += (float(b[j-1])-0.5)*(float(b[k-1])-0.5)*v[i]
                else:
                    continue
    return w

# The transverse part of the square of the total spin operator: \sum_{i=0}^{L-1}\sum_{j=0}^{L-1}1/2(S_{i}^{+}S_{j}^{-}+S_{i}^{-}S_{j}^{+}).

def Transverse(v):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        s = spinBasis[i % subDim] #i = h*dim_sub+s
        h = holeBasis[i // subDim]
        b = bin(s)[2:]
        b = '%0*.0f' % (numSite, int(b))  #the highest bit is filled with a '0' just for convenience while does not have real meaning
        b = b[::-1]
        for j in range(numSite):
            for k in range(numSite):
                if j == k and j != h:
                    w[i] += 0.5*v[i]
                elif j < h and k < h and b[j] != b[k]:
                    ss = Flip(s, j, k)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                elif j > h and k < h and b[j-1] != b[k]:
                    ss = Flip(s, j-1, k)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                elif j < h and k > h and b[j] != b[k-1]:
                    ss = Flip(s, j, k-1)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                elif j > h and k > h and b[j-1] != b[k-1]:
                    ss = Flip(s, j-1, k-1)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                else:
                    continue
    return w

def Region(l, h): # Only valid for the periodic boundary condition.
    s = []
    if l == 0:
        return s
    elif (l % 2):
        halfRight = l//2+1
        halfLeft = l//2
    else:
        halfRight = l//2
        halfLeft = l//2
    for i in range(1, halfRight+1):
        s.append((h+i) % numSite)
    for i in range(1, halfLeft+1):
        s.append((h+numSite-i) % numSite)
    return s

def ProjectedHole(v, hh):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        h = holeBasis[i // subDim]
        if h != hh:
            continue
        else:
            w[i] = v[i]
    return w
 
# To compute the wavefunction projected to a specific hole position hh with length l around hh.

def ProjectedSzSquare(v, hh, l):
    w = np.zeros(dim, dtype=complex)
    reg = Region(l, hh)
#    print(hh, reg)
    for i in range(dim):
        h = holeBasis[i // subDim]
        if h != hh:
            continue
        s = spinBasis[i % subDim] # i = h*subDim+(i % sumDim)
        b = bin(s)[2:]
        b = '%0*.0f' % (numSite, int(b))  # The highest bit is filled with a '0' just for convenience while does not have real meaning
        b = b[::-1]
        for j in reg:
            for k in reg:
                if j < h and k < h:
                    w[i] += (float(b[j])-0.5)*(float(b[k])-0.5)*v[i]
                elif j < h and k > h:
                    w[i] += (float(b[j])-0.5)*(float(b[k-1])-0.5)*v[i]
                elif j > h and k < h:
                    w[i] += (float(b[j-1])-0.5)*(float(b[k])-0.5)*v[i]
                elif j > h and k > h:
                    w[i] += (float(b[j-1])-0.5)*(float(b[k-1])-0.5)*v[i]
    return w

def ProjectedTransverse(v, hh, l):
    w = np.zeros(dim, dtype=complex)
    reg = Region(l, hh)
    for i in range(dim):
        h = holeBasis[i // subDim]
        if h != hh:
            continue
        s = spinBasis[i % subDim] 
        b = bin(s)[2:]
        b = '%0*.0f' % (numSite, int(b))
        b = b[::-1]
        for j in reg:
            for k in reg:
                if j == k:
                    w[i] += 0.5*v[i]
                elif j < h and k < h and b[j] != b[k]:
                    ss = Flip(s, j, k)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                elif j > h and k < h and b[j-1] != b[k]:
                    ss = Flip(s, j-1, k)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                elif j < h and k > h and b[j] != b[k-1]:
                    ss = Flip(s, j, k-1)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                elif j > h and k > h and b[j-1] != b[k-1]:
                    ss = Flip(s, j-1, k-1)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*v[i]
                else:
                    continue
    return w

def SampleProjectedSzSquare(i, h):
    temp = np.zeros((2, 2), dtype=complex)
    m = np.zeros(numSite, dtype=complex)
    if np.absolute(ene[i][0]-ene[i][1]) < 0.0000001:
        wf0t = Translation(wf0[i])
        wf1t = Translation(wf1[i])
        temp[0][0] = np.vdot(wf0[i], wf0t) # vdot for complex conjugate...
        temp[0][1] = np.vdot(wf0[i], wf1t)
        temp[1][0] = np.vdot(wf1[i], wf0t)
        temp[1][1] = np.vdot(wf1[i], wf1t)
        w, v = la.eig(temp)
        print(i, np.angle(w, deg=True))
        if w[0].imag > 0:
            choice = 0
        else:
            choice = 1
        wf = v[:, choice][0]*wf0[i]+v[:, choice][1]*wf1[i]  #construct eigenstates with specific momentum
    else:
        wf = wf0[i]
    print('Total Sz square: ', np.vdot(wf, TotalSzSquare(wf)))
    for l in range(numSite):
        m[l] = np.vdot(wf, ProjectedSzSquare(wf, h, l))
    return np.real(m)

def SampleProjectedTransverse(i, h):
    temp = np.zeros((2, 2), dtype=complex)
    m = np.zeros(numSite, dtype=complex)
    if np.absolute(ene[i][0]-ene[i][1]) < 0.0000001:
        wf0t = Translation(wf0[i])
        wf1t = Translation(wf1[i])
        temp[0][0] = np.vdot(wf0[i], wf0t) # vdot for complex conjugate...
        temp[0][1] = np.vdot(wf0[i], wf1t)
        temp[1][0] = np.vdot(wf1[i], wf0t)
        temp[1][1] = np.vdot(wf1[i], wf1t)
        w, v = la.eig(temp)
        print(i, np.angle(w, deg=True))
        if w[0].imag > 0:
            choice = 0
        else:
            choice = 1
        wf = v[:, choice][0]*wf0[i]+v[:, choice][1]*wf1[i]  #construct eigenstates with specific momentum
    else:
        wf = wf0[i]
    print('Total transerse: ', np.vdot(wf, Transverse(wf)))
    for l in range(numSite):
        m[l] = np.vdot(wf, ProjectedTransverse(wf, h, l))
    return np.real(m)

def Sample(paras, i):
    temp = np.zeros((2, 2), dtype=complex)
    m = np.zeros(numSite, dtype=complex)
    if np.absolute(ene[i][0]-ene[i][1]) < 0.0000000001:
        wf0t = Translation(wf0[i])
        wf1t = Translation(wf1[i])
        temp[0][0] = np.vdot(wf0[i], wf0t) # vdot for complex conjugate...
        temp[0][1] = np.vdot(wf0[i], wf1t)
        temp[1][0] = np.vdot(wf1[i], wf0t)
        temp[1][1] = np.vdot(wf1[i], wf1t)
        w, v = la.eig(temp)
        print(i, np.angle(w, deg=True))
        if w[0].imag > 0: # choose the eige
            choice = 0
        else:
            choice = 1
        wf = v[:, choice][0]*wf0[i]+v[:, choice][1]*wf1[i]  #construct eigenstates with specific momentum
        wff = Translation(wf)
        wff = np.conjugate(w[choice])*wff
        print('Test trsanlation:', wf-wff)
    else:
        wf = wf0[i]
    tot_sz = 0.0
    tot_trans = 0.0
    for l in range(numSite):
        m[l] = np.vdot(wf, ProjectedHole(wf, l))
        tot_sz += np.vdot(wf, ProjectedSzSquare(wf, l, numSite-1))
        tot_trans += np.vdot(wf, ProjectedTransverse(wf, l, numSite-1))
        print(np.vdot(wf, ProjectedSzSquare(wf, l, numSite-1)), np.vdot(wf, ProjectedTransverse(wf, l, numSite-1)))
    print('Total: ', tot_sz, tot_trans)
    return m

def TransverseAlongJ(wf0, wf1):
    temp = np.zeros((2, 2), dtype=complex)
    m = np.zeros(numSam)
    for i in range(numSam):
        if np.absolute(ene[i][0]-ene[i][1]) < 0.0000000001:
            wf0t = Translation(wf0[i])
            wf1t = Translation(wf1[i])
            temp[0][0] = np.vdot(wf0[i], wf0t) # vdot for complex conjugate...
            temp[0][1] = np.vdot(wf0[i], wf1t)
            temp[1][0] = np.vdot(wf1[i], wf0t)
            temp[1][1] = np.vdot(wf1[i], wf1t)
            w, v = la.eig(temp)
            print(i, np.angle(w, deg=True))
            if w[0].imag > 0: # choose the eige
                choice = 0
            else:
                choice = 1
                wf = v[:, choice][0]*wf0[i]+v[:, choice][1]*wf1[i]  #construct eigenstates with specific momentum
        else:
            wf = wf0[i]
        m[i] = np.vdot(wf, Transverse(wf))
        print(i, m[i])
    return m

def NSz(v, h, j):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        hh = holeBasis[i // subDim]
        if h != hh:
            continue
        s = spinBasis[i % subDim] #i = h*dim_sub+s
        b = bin(s)[2:]
        b = '%0*.0f' % (numSite, int(b))  #the highest bit is filled with a '0' just for convenience while does not have real meaning
        b = b[::-1]
        if j < hh:
            w[i] += (float(b[j])-0.5)*v[i]
        elif j > hh:
            w[i] += (float(b[j-1])-0.5)*v[i]
        else:
            continue
    return w
 
def SampleNSz(i, h):
    temp = np.zeros((2, 2), dtype=complex)
    m = np.zeros(numSite, dtype=complex)
    if np.absolute(ene[i][0]-ene[i][1]) < 1e-10:
        wf0t = Translation(wf0[i])
        wf1t = Translation(wf1[i])
        temp[0][0] = np.vdot(wf0[i], wf0t) # vdot for complex conjugate...
        temp[0][1] = np.vdot(wf0[i], wf1t)
        temp[1][0] = np.vdot(wf1[i], wf0t)
        temp[1][1] = np.vdot(wf1[i], wf1t)
        w, v = la.eig(temp)
        print(i, np.angle(w, deg=True))
        if w[0].imag > 0:
            choice = 0
        else:
            choice = 1
        wf = v[:, choice][0]*wf0[i]+v[:, choice][1]*wf1[i] 
    else:
        wf = wf0[i]
    for l in range(numSite):
        m[l] = np.vdot(wf, NSz(wf, h, l))
    return np.real(m)

hole = 5

slice1 = 10
slice2 = 350
paras = (numEval, numSite, numSam, 'PBC', 'false', flux)
print(paras, slice)
ene = LoadEigVal(paras)
wf0 = LoadEigVec(paras, 0)
wf1 = LoadEigVec(paras, 1)
nz1 = SampleNSz(slice1, hole)
nz2 = SampleNSz(slice2, hole)
pz1 = SampleProjectedSzSquare(slice1, hole)
pz2 = SampleProjectedSzSquare(slice2, hole)

paras = (numEval, numSite, numSam, 'PBC', 'true', flux)
print(paras, slice)
ene = LoadEigVal(paras)
wf0 = LoadEigVec(paras, 0)
wf1 = LoadEigVec(paras, 1)
nz1_sigma = SampleNSz(slice1, hole)
nz2_sigma = SampleNSz(slice2, hole)
pz1_sigma = SampleProjectedSzSquare(slice1, hole)
pz2_sigma = SampleProjectedSzSquare(slice2, hole)

fig = plt.figure(figsize=(8, 8))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
x = np.arange(numSite)

yLimLow = -0.02
yLimUp = 0.04

ax1 = plt.subplot(221)
ax1.set_title('(a) $J/t=1.0$', fontsize=fs, x=0.25, y=0.85)
plt.plot(x, nz1, '-o', color='red', label='$t$-$J$')
plt.plot(x, nz1_sigma, '-->', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.axvline(x=5, ymin=0.0, ymax=1.0, color='black', linestyle='-.')
plt.legend(loc='upper right', frameon=False, fontsize=fs)
plt.ylim([yLimLow, yLimUp])
plt.setp(ax1.get_xticklabels(), fontsize=fs)
plt.setp(ax1.get_yticklabels(), fontsize=fs)
plt.yticks(np.arange(yLimLow, yLimUp, 0.01))
plt.xlabel('$j$', fontsize=fs)
plt.ylabel('$C^{z}(h, l)$', fontsize=fs)
# ax1.text(0.1, 0.75, '', fontsize=fs, )

ax2 = plt.subplot(222, sharey=ax1)
ax2.set_title('(b) $J/t=35.0$', fontsize=fs, x=0.25, y=0.85)
plt.plot(x, nz2, '-o', color='red', label='$t$-$J$')
plt.plot(x, nz2_sigma, '-->', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.axvline(x=5, ymin=0.0, ymax=1.0, color='black', linestyle='-.')
plt.legend(loc='upper right', frameon=False, fontsize=fs)
plt.ylim([yLimLow, yLimUp])
plt.yticks(np.arange(yLimLow, yLimUp, 0.01))
plt.setp(ax2.get_xticklabels(), fontsize=fs)
plt.setp(ax2.get_yticklabels(), visible=False, fontsize=fs)
plt.xlabel('$j$', fontsize=fs)
# ax3.text(0.1, 0.75, '', fontsize=fs, transform=ax3.transAxes)

yLimUp = 0.12
yLimLow = -0.06

ax3 = plt.subplot(223, sharex=ax1)
ax3.set_title('(c) $J/t=1.0$', fontsize=fs, x=0.25, y=0.85)
plt.plot(x, pz1, '-o', color='red', label='$t$-$J$') 
plt.plot(x, pz1_sigma, '-->', color='blue', label='$\sigma\cdot{t}$-$J$') 
plt.axhline(y=0.0, xmin=0.0, xmax=1.0, linestyle='-.', color='black')
# plt.plot(x, mark, 'b--')
plt.legend(loc='best', frameon=False, fontsize=fs)
plt.ylim([yLimLow, yLimUp])
plt.yticks(np.arange(yLimLow, yLimUp, 0.03)) # Change the ticks' frequency.
plt.setp(ax3.get_xticklabels(), fontsize=fs)
plt.setp(ax3.get_yticklabels(), fontsize=fs)
plt.xlabel('$l$', fontsize=fs)
plt.ylabel('$R^{z}(h, l)$', fontsize=fs)
# ax1.text(0.2, 0.1, '', fontsize=fs, transform=ax1.transAxes)

ax4 = plt.subplot(224, sharey=ax3, sharex=ax2)
ax4.set_title('(d) $J/t=35.0$', fontsize=fs, x=0.25, y=0.85)
plt.plot(x, pz2, '-o', color='red', label='$t$-$J$') 
plt.plot(x, pz2_sigma, '-->', color='blue', label='$\sigma\cdot{t}$-$J$') 
plt.axhline(y=0.0, xmin=0.0, xmax=1.0, linestyle='-.', color='black')
plt.legend(loc='best', frameon=False, fontsize=fs)
plt.ylim([yLimLow, yLimUp])
plt.yticks(np.arange(yLimLow, yLimUp, 0.03)) # Change the ticks' frequency.
plt.setp(ax4.get_xticklabels(), fontsize=fs)
plt.setp(ax4.get_yticklabels(), visible=False, fontsize=fs)
plt.xlabel('$l$', fontsize=fs)
# ax2.text(0.2, 0.1, '', fontsize=fs, transform=ax3.transAxes)

image = 'spinz_slice_len_%s_slice_%s_flux_%s.pdf'
imageParas = (numSite, slice, flux)
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
plt.show()
