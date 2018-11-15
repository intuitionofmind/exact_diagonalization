#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from auxi import numSite, numEval, cut, fs

deltaT = 0.001
numTime = 400
inter = 4
numPlot = int(numTime/inter)

plotInter = 2

level = 240

# dataDir = '/Users/wayne/Downloads/data/'
dataDir = './'

meeFile = 'time_evolution_mEE%s_level%s_len%s_J%s_%s_sigma%s.npy'
seeFile = 'time_evolution_sEE%s_level%s_len%s_J%s_%s_sigma%s.npy'

fig = plt.figure(figsize=(10, 10))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
x = np.linspace(0, deltaT*numTime, num=numPlot)


J = 0.3
yt = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

npyParas = (numTime, level,  numSite, J, 'PBC', 'False')
mee = np.load(os.path.join(dataDir, meeFile) % npyParas)
see = np.load(os.path.join(dataDir, seeFile) % npyParas)

ax1 = plt.subplot(221)
ax1.set_title('(a) $t$-$J$, $J=0.3$', fontsize=fs, x=0.2, y=0.1)
plt.plot(x[::plotInter], mee[::plotInter], '.-', color='red', label='mEE')
plt.plot(x[::plotInter], see[::plotInter], 'x-', color='blue', label='bEE')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.setp(ax1.get_xticklabels(), visible=False, fontsize=fs)
ax1.set_xlim([0, 0.4])
ax1.set_ylim([0.0, 3.5])
plt.yticks(yt, fontsize=fs)
print(np.mean(mee[int(numPlot/2):]), np.std(mee[int(numPlot/2):])/(np.sqrt(len(mee)/2)))
print(np.mean(see[int(numPlot/2):]), np.std(see[int(numPlot/2):])/(np.sqrt(len(see)/2)))

npyParas = (numTime, level, numSite, J, 'PBC', 'True')
mee = np.load(os.path.join(dataDir, meeFile) % npyParas)
see = np.load(os.path.join(dataDir, seeFile) % npyParas)

ax2 = plt.subplot(222, sharey=ax1)
ax2.set_title('(b) $\sigma\cdot{t}$-$J$, $J=0.3$', fontsize=fs, x=0.25, y=0.1)
plt.plot(x[::plotInter], mee[::plotInter], '.-', color='red', label='mEE')
plt.plot(x[::plotInter], see[::plotInter], 'x-', color='blue', label='bEE')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.setp(ax2.get_xticklabels(), visible=False, fontsize=fs)
plt.setp(ax2.get_yticklabels(), visible=False, fontsize=fs)
ax2.set_xlim([0, 0.4])
ax2.set_ylim([0, 3.5])
print(np.mean(mee[int(numPlot/2):]), np.std(mee[int(numPlot/2):])/(np.sqrt(len(mee)/2)))
print(np.mean(see[int(numPlot/2):]), np.std(see[int(numPlot/2):])/(np.sqrt(len(see)/2)))

J = 40.3
xt = [0.0, 0.1, 0.2, 0.3, 0.4]

npyParas = (numTime, level, numSite, J, 'PBC', 'False')
mee = np.load(os.path.join(dataDir, meeFile) % npyParas)
see = np.load(os.path.join(dataDir, seeFile) % npyParas)

ax3 = plt.subplot(223, sharex=ax1)
ax3.set_title('(c) $t$-$J$, $J=40.3$', fontsize=fs, x=0.2, y=0.1)
plt.plot(x[::plotInter], mee[::plotInter], '.-', color='red', label='mEE')
plt.plot(x[::plotInter], see[::plotInter], 'x-', color='blue', label='bEE')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.xlabel('$Tt$', fontsize=fs)
plt.xticks(xt, fontsize=fs)
plt.yticks(yt, fontsize=fs)
ax3.set_xlim([0, 0.4])
ax3.set_ylim([0, 3.5])
print(np.mean(mee[int(numPlot/2):]), np.std(mee[int(numPlot/2):])/(np.sqrt(len(mee)/2)))
print(np.mean(see[int(numPlot/2):]), np.std(see[int(numPlot/2):])/(np.sqrt(len(see)/2)))

npyParas = (numTime, level, numSite, J, 'PBC', 'True')
mee = np.load(os.path.join(dataDir, meeFile) % npyParas)
see = np.load(os.path.join(dataDir, seeFile) % npyParas)

ax4 = plt.subplot(224, sharey=ax3)
ax4.set_title('(d) $\sigma\cdot{t}$-$J$, $J=40.3$', fontsize=fs, x=0.25, y=0.1)
plt.plot(x[::plotInter], mee[::plotInter], '.-', color='red', label='mEE')
plt.plot(x[::plotInter], see[::plotInter], 'x-', color='blue', label='bEE')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.setp(ax4.get_yticklabels(), visible=False, fontsize=fs)
plt.xlabel('$Tt$', fontsize=fs)
plt.xticks(xt, fontsize=fs)
ax4.set_xlim([0, 0.4])
ax4.set_ylim([0, 3.5])
print(np.mean(mee[int(numPlot/2):]), np.std(mee[int(numPlot/2):])/(np.sqrt(len(mee)/2)))
print(np.mean(see[int(numPlot/2):]), np.std(see[int(numPlot/2):])/(np.sqrt(len(see)/2)))

image = 'time_evolution_EE_len%s_%s.pdf'
imageParas = (numSite, 'PBC')
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
plt.show()
