#!/usr/bin/python3

import os
import math
import random
import numpy as np
import matplotlib as mpl
from scipy import linalg as la
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from tjchain import LoadEigVal
from tjchain import cut, numSite, numEval, numSam

fs = 24 # font size 

dataDir = '/Users/wayne/Downloads/data/'

meeFile = 'mutualEE%s_len%s_J%s_%s_sigma%s.npy'
mesFile = 'mutualES%s_len%s_J%s_%s_sigma%s.npy'
seeFile = 'spatialEE%s_len%s_cut%s_J%s_%s_sigma%s.npy'
sesFile = 'spatialES%s_len%s_cut%s_J%s_%s_sigma%s.npy'

fig = plt.figure(figsize=(10, 10))
mpl.rcParams['axes.linewidth'] = 2.0
plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

yt = [1.0, 2.0, 3.0, 4.0]

eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'
l = 0
J = 0.3

mEEnpyParasF = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParasF = (numEval, numSite, cut, J, 'PBC', 'False')
mEEnpyParasT = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParasT = (numEval, numSite, cut, J, 'PBC', 'True')

ax1 = plt.subplot(221)
paras = (numEval, numSite, numSam, 'PBC', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
print('ground state energy', ene[l][0])
x = (ene[l]-ene[l][0])/numSite
print(ene[l])

meeF = np.load(os.path.join(dataDir, meeFile) % mEEnpyParasF)
meeT = np.load(os.path.join(dataDir, meeFile) % mEEnpyParasT)
# plt.plot(x[::inter], ee[::inter], '-o', color='red', label='$t$-$J$')
plt.plot(x, meeF, '.', color='red', label='$t$-$J$')
plt.plot(x, meeT, 'x', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.legend(loc='upper right', frameon=False, fontsize=fs)

ax1.set_xlim(0.0, 0.5)
ax1.set_ylim(0.0, 5.0)

xstart, xend = ax1.get_xlim()
ax1.xaxis.set_ticks(np.arange(xstart, xend, 0.1))
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ystart, yend = ax1.get_ylim()
ax1.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax1.tick_params(axis='both', labelsize=fs, direction='in')

plt.ylabel('EE', fontsize=fs)
ax1.set_title('(a) mEE, $J/t=0.3$', fontsize=fs, x=0.35, y=0.9)

ax2 = plt.subplot(222, sharey=ax1)
paras = (numEval, numSite, numSam, 'PBC', 'True')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
x = (ene[l]-ene[l][0])/numSite

# mEEnpyParas = (numEval, numSite, J, 'PBC', 'True')
# sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
# mee = np.load(os.path.join(dataDir, meeFile) % mEEnpyParas)
seeF = np.load(os.path.join(dataDir, seeFile) % sEEnpyParasF)
seeT = np.load(os.path.join(dataDir, seeFile) % sEEnpyParasT)
plt.plot(x, seeF, '.', color='red', label='$t$-$J$')
plt.plot(x, seeT, 'x', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.legend(loc='upper right', frameon=False, fontsize=fs)

ax2.set_xlim(0.0, 0.5)
ax2.set_ylim(0.0, 5.0)

xstart, xend = ax2.get_xlim()
ax2.xaxis.set_ticks(np.arange(xstart, xend, 0.1))
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ystart, yend = ax2.get_ylim()
ax2.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax2.tick_params(axis='both', labelsize=fs, direction='in')
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.set_title('(b) bEE, $J/t=0.3$', fontsize=fs, x=0.35, y=0.9)

l = 4
J = 40.3

mEEnpyParasF = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParasF = (numEval, numSite, cut, J, 'PBC', 'False')
mEEnpyParasT = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParasT = (numEval, numSite, cut, J, 'PBC', 'True')

ax3 = plt.subplot(223)
paras = (numEval, numSite, numSam, 'PBC', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
x = (ene[l]-ene[l][0])/numSite

meeF = np.load(os.path.join(dataDir, meeFile) % mEEnpyParasF)
meeT = np.load(os.path.join(dataDir, meeFile) % mEEnpyParasT)
# plt.plot(x[::inter], ee[::inter], '-o', color='red', label='$t$-$J$')
plt.plot(x, meeF, '.', color='red', label='$t$-$J$')
plt.plot(x, meeT, 'x', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.legend(loc='upper right', frameon=False, fontsize=fs)

ax3.set_xlim(0.0, 20.0)
ax3.set_ylim(0.0, 5.0)
xstart, xend = ax3.get_xlim()
ax3.xaxis.set_ticks(np.arange(xstart, xend, 5.0))
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ystart, yend = ax3.get_ylim()
ax3.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax3.tick_params(axis='both', labelsize=fs, direction='in')
plt.xlabel('$E/L$', fontsize=fs)
plt.ylabel('EE', fontsize=fs)
ax3.set_title('(c) mEE, $J/t=40.3$', fontsize=fs, x=0.35, y=0.9)

ax4 = plt.subplot(224, sharey=ax3)
ax4.set_xlim(0.0, 20.0)
ax4.set_ylim(0.0, 5.0)

paras = (numEval, numSite, numSam, 'PBC', 'True')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile % paras))
x = (ene[l]-ene[l][0])/numSite

seeF = np.load(os.path.join(dataDir, seeFile) % sEEnpyParasF)
seeT = np.load(os.path.join(dataDir, seeFile) % sEEnpyParasT)
plt.plot(x, seeF, '.', color='red', label='$t$-$J$')
plt.plot(x, seeT, 'x', color='blue', label='$\sigma\cdot{t}$-$J$')
plt.legend(loc='upper right', frameon=False, fontsize=fs)

xstart, xend = ax4.get_xlim()
ax4.xaxis.set_ticks(np.arange(xstart, xend, 5.0))
ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ystart, yend = ax4.get_ylim()
ax4.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax4.tick_params(axis='both', labelsize=fs, direction='in')
plt.setp(ax4.get_yticklabels(), visible=False)
plt.xlabel('$E/L$', fontsize=fs)
ax4.set_title('(d) bEE, $J/t=40.3$', fontsize=fs, x=0.35, y=0.9)

image = 'excited_states_EE_2.pdf'
# imageParas = (numSite, cut, 'PBC')
fig.tight_layout()
plt.savefig(image, format='pdf')
# plt.show()
