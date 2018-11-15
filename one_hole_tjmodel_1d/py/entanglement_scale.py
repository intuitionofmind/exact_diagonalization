#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from auxi import LoadEigVal
from auxi import Hamiltonian
from auxi import cut, numSite, numEval, numSam

fs = 16 # font size 

dataDir = '/Users/wayne/Downloads/data/'

eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'

meeFile = 'mutualEE%s_len%s_J%s_%s_sigma%s.npy'
mesFile = 'mutualES%s_len%s_J%s_%s_sigma%s.npy'
seeFile = 'spatialEE%s_len%s_cut%s_J%s_%s_sigma%s.npy'
sesFile = 'spatialES%s_len%s_cut%s_J%s_%s_sigma%s.npy'

def MeasureSEE(paras, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    see = np.load(os.path.join(dataDir, seeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[l][i])/paras[1]
        if  e > (targetE-deltaE) and e < (targetE+deltaE):
            temp = np.append(temp, see[i])
    print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))

def MeasureMEE(paras, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    mee = np.load(os.path.join(dataDir, meeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[l][i])/paras[1]
        if e > (targetE-deltaE) and e < (targetE+deltaE):
            temp = np.append(temp, mee[i])
    print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))

# print(MeasureSEE(sEEnpyParas, targetE, deltaE), MeasureMEE(mEEnpyParas, targetE, deltaE))

fig = plt.figure(figsize=(10, 10))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
x = [10, 12, 14, 16]
xLim = [9.8, 16.2]

num = 4
meeScale = np.zeros(num)
meeScaleErr = np.zeros(num)
seeScale = np.zeros(num)
seeScaleErr = np.zeros(num)

level = 240 # excited eigenstate number 

cut = 4

l = 0
J = 0.3

deltaE = 0.004 # Resolution of the energy density

yt = [1.0, 1.5, 2.0, 2.5]
yLim = [0.5, 3.0]

for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'False')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
    # print(eneSpectrum[l][level]-eneSpectrum[l][0])
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'False')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale[i] = res[0]
    meeScaleErr[i] = res[1]
    sEEnpyParas = (numEval, size, cut, J, 'PBC', 'False')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale[i] = res[0]
    seeScaleErr[i] = res[1]

ax1 = plt.subplot(221)
ax1.set_title('(a) $t$-$J$, $J=0.3$', fontsize=fs, x=0.2, y=0.9)
# p1 = np.polyfit(xx, seeScale, 1)
# ax1.plot(xx, xx*p1[0]+p1[1], linestyle='--', color='blue')
ax1.errorbar(x, meeScale, yerr=meeScaleErr, linestyle='--', marker='>', linewidth=1.0, color='blue', label='mEE')
# ax1.errorbar(x, seeScale, yerr=seeScaleErr, linestyle='--', marker='.', linewidth=2.0, color='blue', label='bEE')
# plt.legend(loc='upper left', frameon=False, fontsize=fs)
ax1.set_xlim(xLim)
ax1.set_ylim(yLim)
plt.ylabel('mEE', fontsize=fs)
plt.xticks(x, fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.setp(ax1.get_xticklabels(), visible=False)

for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'False')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
    # print(eneSpectrum[l][level]-eneSpectrum[l][0])
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'True')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale[i] = res[0]
    meeScaleErr[i] = res[1]
    sEEnpyParas = (numEval, size, cut, J, 'PBC', 'True')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale[i] = res[0]
    seeScaleErr[i] = res[1]

ax2 = plt.subplot(222, sharey=ax1)
ax2.set_title('(b) $\sigma\cdot{t}$-$J$, $J=0.3$', fontsize=fs, x=0.25, y=0.9)
# p2 = np.polyfit(xx, seeScale, 1)
# ax2.plot(xx, xx*p2[0]+p2[1], linestyle='--', color='blue')
ax2.errorbar(x, meeScale, yerr=meeScaleErr, linestyle='--', marker='>', linewidth=1.0, color='blue', label='mEE')
# ax2.errorbar(x, seeScale, yerr=seeScaleErr, linestyle=' ', marker='.', linewidth=2.0, color='blue', label='bEE')
# plt.legend(loc='upper left', frameon=False, fontsize=fs)
ax2.set_xlim(xLim)
ax2.set_ylim(yLim)
plt.xticks(x, fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

l=4
J=40.3
deltaE = 2.0 # Resolution of the energy density

yLim = [1.6, 2.2]
yt = [1.7, 1.8, 1.9, 2.0, 2.1]

for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'False')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
    # print(eneSpectrum[l][level]-eneSpectrum[l][0])
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'False')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale[i] = res[0]
    meeScaleErr[i] = res[1]
    sEEnpyParas = (numEval, size, cut, J, 'PBC', 'False')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale[i] = res[0]
    seeScaleErr[i] = res[1]

ax3 = plt.subplot(223, sharex=ax1)
ax3.set_title('(c) $t$-$J$, $J=40.3$', fontsize=fs, x=0.2, y=0.9)
# p3 = np.polyfit(xx, seeScale, 1)
# ax3.plot(xx, xx*p3[0]+p3[1], linestyle='--', color='blue')
ax3.errorbar(x, meeScale, yerr=meeScaleErr, linestyle='--', marker='>', linewidth=1.0, color='blue', label='mEE')
# ax3.errorbar(x, seeScale, yerr=seeScaleErr, linestyle=' ', marker='.', linewidth=2.0, color='blue', label='bEE')
# plt.legend(loc='upper left', frameon=False, fontsize=fs)
ax3.set_xlim(xLim)
ax3.set_ylim(yLim)
plt.xticks(x, fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.xlabel('L', fontsize=fs)
plt.ylabel('mEE', fontsize=fs)

for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'True')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
    # print(eneSpectrum[l][level]-eneSpectrum[l][0])
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'True')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale[i] = res[0]
    meeScaleErr[i] = res[1]
    sEEnpyParas = (numEval, size, cut, J, 'PBC', 'True')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale[i] = res[0]
    seeScaleErr[i] = res[1]

ax4 = plt.subplot(224, sharex=ax2, sharey=ax3)
ax4.set_title('(d) $\sigma\cdot{t}$-$J$, $J=40.3$', fontsize=fs, x=0.25, y=0.9)
# p4 = np.polyfit(xx, seeScale, 1)
# ax4.plot(xx, xx*p4[0]+p4[1], linestyle='--', color='blue')
ax4.errorbar(x, meeScale, yerr=meeScaleErr, linestyle='--', marker='>', linewidth=1.0, color='blue', label='mEE')
# ax4.errorbar(x, seeScale, yerr=seeScaleErr, linestyle=' ', marker='.', linewidth=2.0, color='blue', label='bEE')
# plt.legend(loc='upper left', frameon=False, fontsize=fs)
ax4.set_xlim(xLim)
ax4.set_ylim(yLim)
plt.xticks(x, fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.xlabel('L', fontsize=fs)
plt.setp(ax4.get_yticklabels(), visible=False)

image = 'mutual_EE_scale_%s.pdf'
imageParas = ('PBC')
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
plt.show()
