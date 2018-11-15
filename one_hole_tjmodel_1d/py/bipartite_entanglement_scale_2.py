#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib as mpl
from auxi_1d import LoadEigVal
# from auxi_1d import Hamiltonian
from auxi_1d import cut, numSite, numEval, numSam

fs = 20 # font size 

dataDir = '/Users/wayne/Downloads/data/'
eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'
seeFile = 'spatialEE%s_len%s_cut%s_J%s_%s_sigma%s.npy'
sesFile = 'spatialES%s_len%s_cut%s_J%s_%s_sigma%s.npy'

def MeasureSEE(paras, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    see = np.load(os.path.join(dataDir, seeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[l][i])/numSite
        if  e > (targetE-deltaE) and e < (targetE+deltaE):
            temp = np.append(temp, see[i])
    print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))

def MeasureMEE(paras, targetE, deltaE):
    temp = np.zeros(0, dtype=float)
    mee = np.load(os.path.join(dataDir, meeFile) % paras)
    for i in range(numEval):
        e = (eneSpectrum[l][i])/numSite
        if e > (targetE-deltaE) and e < (targetE+deltaE):
            temp = np.append(temp, mee[i])
    print(len(temp))
    return (np.mean(temp), np.std(temp)/np.sqrt(len(temp)))

# print(MeasureSEE(sEEnpyParas, targetE, deltaE), MeasureMEE(mEEnpyParas, targetE, deltaE))

fig = plt.figure(figsize=(8, 4))
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
x = [2, 3, 4, 5]
xx = np.log(x)

num = 4
seeScale0 = np.zeros(num)
seeScaleErr0 = np.zeros(num)
seeScale1 = np.zeros(num)
seeScaleErr1 = np.zeros(num)

level = 299 # excited eigenstate number 

l = 0
J = 0.3

paras = (numEval, numSite, numSam, 'PBC', 'False')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
targetE = (eneSpectrum[l][level])/numSite
deltaE = 0.01 # Resolution of the energy density
for i in range(num):
    cut = i+2
    sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale0[i] = res[0]
    seeScaleErr0[i] = res[1]

paras = (numEval, numSite, numSam, 'PBC', 'True')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
print(eneSpectrum[l][level]-eneSpectrum[l][0])
targetE = (eneSpectrum[l][level])/numSite
deltaE = 0.01 # Resolution of the energy density
for i in range(num):
    cut = i+2
    sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale1[i] = res[0]
    seeScaleErr1[i] = res[1]

xLim = [1.9, 5.2]
yt = [2.0, 2.5, 3.0, 3.5]
yLim = [1.6, 3.3]

ax1 = plt.subplot(121)
ax1.set_title('(a) $J/t=0.3$', fontsize=fs, x=0.3, y=0.85)
p0 = np.polyfit(xx, seeScale0, 1)
p1 = np.polyfit(xx, seeScale1, 1)
ax1.errorbar(xx, seeScale0, yerr=seeScaleErr0, linestyle=' ', marker='x', linewidth=2.0, color='blue', label='$t$-$J$')
ax1.plot(xx, xx*p0[0]+p0[1], linestyle='--', color='blue')
ax1.errorbar(xx, seeScale1, yerr=seeScaleErr1, linestyle=' ', marker='>', linewidth=2.0, color='red', label='$\sigma\cdot{t}$-$J$')
ax1.plot(xx, xx*p1[0]+p1[1], linestyle='--', color='red')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.xticks(xx, x, fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.xlabel('Bipartite cut', fontsize=fs)
plt.ylabel('bEE', fontsize=fs)
ax1.set_xlim(np.log(xLim))
ax1.set_ylim(yLim)
ax1.tick_params(axis='both', labelsize=fs, direction='in')
# plt.setp(ax1.get_xticklabels(), visible=False)


l=4
J=40.3

paras = (numEval, numSite, numSam, 'PBC', 'False')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
# print(eneSpectrum[l][level]-eneSpectrum[l][0])
targetE = (eneSpectrum[l][level])/numSite
deltaE = 0.2 # Resolution of the energy density
for i in range(num):
    cut = i+2
    sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale0[i] = res[0]
    seeScaleErr0[i] = res[1]

paras = (numEval, numSite, numSam, 'PBC', 'True')
eneSpectrum = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
print(eneSpectrum[l][level]-eneSpectrum[l][0])
targetE = (eneSpectrum[l][level])/numSite
deltaE = 0.2 # Resolution of the energy density
for i in range(num):
    cut = i+2
    sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
    res = MeasureSEE(sEEnpyParas, targetE, deltaE)
    seeScale1[i] = res[0]
    seeScaleErr1[i] = res[1]


xLim = [1.9, 5.2]
yt = [2.0, 2.5, 3.0, 3.5]
yLim = [1.6, 3.5]

ax2 = plt.subplot(122, sharey=ax1)
ax2.set_title('(b) $J/t=40.3$', fontsize=fs, x=0.3, y=0.85)
p0 = np.polyfit(xx, seeScale0, 1)
p1 = np.polyfit(xx, seeScale1, 1)
ax2.errorbar(xx, seeScale0, yerr=seeScaleErr0, linestyle=' ', marker='x', linewidth=2.0, color='blue', label='$t$-$J$')
ax2.plot(xx, xx*p0[0]+p0[1], linestyle='--', color='blue')
ax2.errorbar(xx, seeScale1, yerr=seeScaleErr1, linestyle=' ', marker='>', linewidth=2.0, color='red', label='$\sigma\cdot{t}$-$J$')
ax2.plot(xx, xx*p1[0]+p1[1], linestyle='--', color='red')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.xticks(xx, x, fontsize=fs)
# plt.yticks(yt, fontsize=fs)
plt.xlabel('Bipartite cut', fontsize=fs)
# plt.ylabel('bEE', fontsize=fs)
ax2.set_xlim(np.log(xLim))
# ax2.set_ylim(yLim)
ax2.tick_params(axis='both', labelsize=fs, direction='in')
# plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

image = 'bipartite_EE_scale_len%s_%s.pdf'
imageParas = (numSite, 'PBC')
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
# plt.show()
