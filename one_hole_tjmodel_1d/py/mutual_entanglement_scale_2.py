#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from auxi_1d import LoadEigVal
from auxi_1d import cut, numSite, numEval, numSam

fs = 20 # font size 

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

fig = plt.figure(figsize=(8, 4))
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
x = [10, 12, 14, 16]
xLim = [9.8, 16.2]

num = 4
meeScale0 = np.zeros(num)
meeScaleErr0 = np.zeros(num)

meeScale1 = np.zeros(num)
meeScaleErr1 = np.zeros(num)

level = 299 # excited eigenstate number 

cut = 4

yt = [1.0, 1.5, 2.0, 2.5]
yLim = [0.5, 3.0]

l = 0
J = 0.3
deltaE = 0.004 # Resolution of the energy density
for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'False')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'False')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale0[i] = res[0]
    meeScaleErr0[i] = res[1]

for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'True')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'True')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale1[i] = res[0]
    meeScaleErr1[i] = res[1]

ax1 = plt.subplot(121)
ax1.set_title('(a)$J/t=0.3$', fontsize=fs, x=0.3, y=0.85)
ax1.errorbar(x, meeScale0, yerr=meeScaleErr0, linestyle='--', marker='>', linewidth=1.0, color='blue', label='$t$-$J$')
ax1.errorbar(x, meeScale1, yerr=meeScaleErr1, linestyle='--', marker='x', linewidth=1.0, color='red', label='$\sigma\cdot{t}$-$J$')
plt.legend(loc='best', frameon=False, fontsize=fs)
ax1.set_xlim(xLim)
ax1.set_ylim(yLim)
plt.xlabel('$L$', fontsize=fs)
plt.ylabel('mEE', fontsize=fs)
plt.xticks(x, fontsize=fs)
plt.yticks(yt, fontsize=fs)
ax1.tick_params(axis='both', labelsize=fs, direction='in')
# plt.setp(ax1.get_xticklabels(), visible=False)

l=4
J=40.3
deltaE = 2.0 # Resolution of the energy density
for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'False')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'False')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale0[i] = res[0]
    meeScaleErr0[i] = res[1]

for i in range(num):
    size = numSite+2*i

    paras = (numEval, size, numSam, 'PBC', 'True')
    eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
    targetE = (eneSpectrum[l][level])/size

    mEEnpyParas = (numEval, size, J, 'PBC', 'True')
    res = MeasureMEE(mEEnpyParas, targetE, deltaE)
    meeScale1[i] = res[0]
    meeScaleErr1[i] = res[1]

yLim = [1.6, 2.2]
yt = [1.7, 1.8, 1.9, 2.0, 2.1]
ax2 = plt.subplot(122, sharey=ax1)
ax2.set_title('(b) $J/t=40.3$', fontsize=fs, x=0.3, y=0.85)
ax2.errorbar(x, meeScale0, yerr=meeScaleErr0, linestyle='--', marker='>', linewidth=1.0, color='blue', label='$t$-$J$')
ax2.errorbar(x, meeScale1, yerr=meeScaleErr1, linestyle='--', marker='x', linewidth=1.0, color='red', label='$\sigma\cdot{t}$-$J$')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
ax2.set_xlim(xLim)
# ax2.set_ylim(yLim)
plt.xticks(x, fontsize=fs)
# plt.yticks(yt, fontsize=fs)
plt.xlabel('$L$', fontsize=fs)
ax2.tick_params(axis='both', labelsize=fs, direction='in')
# plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

l=image = 'mutual_EE_scale_%s.pdf'
imageParas = ('PBC')
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
# plt.show()
