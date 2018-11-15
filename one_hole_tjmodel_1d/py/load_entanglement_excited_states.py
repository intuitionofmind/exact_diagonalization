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

meeFile = 'mutualEE%s_len%s_J%s_%s_sigma%s.npy'
mesFile = 'mutualES%s_len%s_J%s_%s_sigma%s.npy'
seeFile = 'spatialEE%s_len%s_cut%s_J%s_%s_sigma%s.npy'
sesFile = 'spatialES%s_len%s_cut%s_J%s_%s_sigma%s.npy'

fig = plt.figure(figsize=(10, 10))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

yt = [1.0, 2.0, 3.0, 4.0]

eigValsFile = 'eigenvalues%s_size%s_J0.3_10.0_%s_%s_sigma%s.dat'
l = 0
J = 0.3

ax1 = plt.subplot(221)
paras = (numEval, numSite, numSam, 'PBC', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
print('ground state energy', ene[l][0])
x = (ene[l]-ene[l][0])/numSite
print(ene[l])

mEEnpyParas = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
mee = np.load(os.path.join(dataDir, meeFile) % mEEnpyParas)
see = np.load(os.path.join(dataDir, seeFile) % sEEnpyParas)
# plt.plot(x[::inter], ee[::inter], '-o', color='red', label='$t$-$J$')
plt.plot(x, mee, '.', color='red', label='mEE')
plt.plot(x, see, 'x', color='blue', label='bEE')
plt.legend(loc='upper right', frameon=False, fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.ylabel('EE', fontsize=fs)
ax1.set_title('(a) $t$-$J$, $J=0.3$', fontsize=fs, x=0.2, y=0.9)
ax1.set_xlim([0., 0.5])
ax1.set_ylim([0., 4.5])

ax2 = plt.subplot(222, sharey=ax1)
paras = (numEval, numSite, numSam, 'PBC', 'True')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
x = (ene[l]-ene[l][0])/numSite

mEEnpyParas = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
mee = np.load(os.path.join(dataDir, meeFile) % mEEnpyParas)
see = np.load(os.path.join(dataDir, seeFile) % sEEnpyParas)
plt.plot(x, mee, '.', color='red', label='mEE')
plt.plot(x, see, 'x', color='blue', label='bEE')
plt.legend(loc='upper right', frameon=False, fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.set_title('(b) $\sigma\cdot{t}$-$J$, $J=0.3$', fontsize=fs, x=0.25, y=0.9)
ax2.set_xlim([0., 0.5])
ax2.set_ylim([0., 4.5])

l = 4
J = 40.3

ax3 = plt.subplot(223)
paras = (numEval, numSite, numSam, 'PBC', 'False')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
x = (ene[l]-ene[l][0])/numSite

mEEnpyParas = (numEval, numSite, J, 'PBC', 'False')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'False')
mee = np.load(os.path.join(dataDir, meeFile) % mEEnpyParas)
see = np.load(os.path.join(dataDir, seeFile) % sEEnpyParas)
# plt.plot(x[::inter], ee[::inter], '-o', color='red', label='$t$-$J$')
plt.plot(x, mee, '.', color='red', label='mEE')
plt.plot(x, see, 'x', color='blue', label='bEE')
plt.legend(loc='upper right', frameon=False, fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.xlabel('$E/L$', fontsize=fs)
plt.ylabel('EE', fontsize=fs)
ax3.set_title('(c) $t$-$J$, $J=40.3$', fontsize=fs, x=0.2, y=0.9)
ax3.set_xlim([0., 22.0])
ax3.set_ylim([0., 4.5])

ax4 = plt.subplot(224, sharey=ax3)
paras = (numEval, numSite, numSam, 'PBC', 'True')
ene = LoadEigVal(os.path.join(dataDir, eigValsFile), paras)
x = (ene[l]-ene[l][0])/numSite

mEEnpyParas = (numEval, numSite, J, 'PBC', 'True')
sEEnpyParas = (numEval, numSite, cut, J, 'PBC', 'True')
mee = np.load(os.path.join(dataDir, meeFile) % mEEnpyParas)
see = np.load(os.path.join(dataDir, seeFile) % sEEnpyParas)
plt.plot(x, mee, '.', color='red', label='mEE')
plt.plot(x, see, 'x', color='blue', label='bEE')
plt.legend(loc='upper right', frameon=False, fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(yt, fontsize=fs)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.xlabel('$E/L$', fontsize=fs)
ax4.set_title('(d) $\sigma\cdot{t}$-$J$, $J=40.3$', fontsize=fs, x=0.25, y=0.9)
ax4.set_xlim([0., 22.0])
ax4.set_ylim([0., 4.5])

image = 'excited_states_EE_len%s_cut%s_%s.pdf'
imageParas = (numSite, cut, 'PBC')
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
plt.show()
