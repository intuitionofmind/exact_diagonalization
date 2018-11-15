#!/usr/bin/python3

import os
import math
import numpy as np
import matplotlib.pyplot as plt
# from auxi import LoadEigVal
from auxi import *

step = 0.1
numEval = 10
flux = 0.0

inter=10

x = np.linspace(step, step*numSam, num=numSam)

fig = plt.figure(figsize = (8, 6))
plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')

ax1 = plt.subplot(121)
paras = (numEval, numSite, numSam, 'PBC', 'false')
ene = LoadEigVal(paras)
print(ene[::inter, 0])
print(ene[::inter, 1])
plt.plot(x[::inter], ene[::inter, 0]-ene[::inter, 0], '-o', color='red', label='$0$')
plt.plot(x[::inter], ene[::inter, 1]-ene[::inter, 0], '->', color='blue', label='$1$')
plt.plot(x[::inter], ene[::inter, 2]-ene[::inter, 0], '-d', color='green', label='$2$')
plt.plot(x[::inter], ene[::inter, 3]-ene[::inter, 0], '-x', color='yellow', label='$3$')
plt.plot(x[::inter], ene[::inter, 4]-ene[::inter, 0], '-h', color='black', label='$4$')
plt.axvline(x=16.9, ymin=-0.15, ymax=1.0, color='blue', linestyle='--')
plt.axvline(x=28.5, ymin=-0.15, ymax=1.0, color='blue', linestyle='--')
plt.xlabel('$J/t$', fontsize=fs)
plt.ylabel('$E-E_{0}~(t)$', fontsize=fs)
plt.ylim([-0.15, 0.5])
plt.legend(loc='best', frameon=False, fontsize=fs)
#plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), fontsize=fs)
ax1.set_title('(a) $t$-$J$', fontsize=fs, x=0.15, y=0.05)

'''
ax1.annotate('CP1', 
        xy=(16.9, -0.02), xycoords='data',
        xytext=(0, -20), textcoords='offset points',
        size=12,
        bbox=dict(boxstyle='round', facecolor=(0.8, 0.9, 0.9), edgecolor='none',),
        arrowprops=dict(arrowstyle='wedge, tail_width=1.0', 
            facecolor=(0.8, 0.9, 0.9), edgecolor='none',
            patchA=None,
            patchB=None,
            relpos=(0.3, 0.5),
            )
        )

ax1.annotate('CP2', 
        xy=(28.5, -0.02), xycoords='data',
        xytext=(0, -20), textcoords='offset points',
        size=12,
        bbox=dict(boxstyle='round', facecolor=(0.8, 0.9, 0.9), edgecolor='none',),
        arrowprops=dict(arrowstyle='wedge, tail_width=1.0', 
            facecolor=(0.8, 0.9, 0.9), edgecolor='none',
            patchA=None,
            patchB=None,
            relpos=(0.3, 0.5),
            )
        )
'''

x = np.linspace(step, step*numSam, num=numSam)

ax2 = plt.subplot(122)
paras = (numEval, numSite, numSam, 'PBC', 'true')
ene = LoadEigVal(paras)
print(ene[0])
print(ene[1])
print(paras)
plt.plot(x[::inter], ene[::inter, 0]-ene[::inter, 0], '-o', color='red', label='$0$')
plt.plot(x[::inter], ene[::inter, 1]-ene[::inter, 0], '->', color='blue', label='$1$')
plt.plot(x[::inter], ene[::inter, 2]-ene[::inter, 0], '-d', color='green', label='$2$')
plt.plot(x[::inter], ene[::inter, 3]-ene[::inter, 0], '-x', color='yellow', label='$3$')
plt.plot(x[::inter], ene[::inter, 4]-ene[::inter, 0], '-h', color='black', label='$4$')
plt.xlabel('$J/t$', fontsize=fs)
plt.ylim([-0.2, 1.2])
plt.legend(loc='best', frameon=False, fontsize=fs)
# plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.set_title('(b) $\sigma\cdot{t}$-$J$', fontsize=fs, x=0.2, y=0.05)

'''
x = np.linspace(step, step*numSam, num=numSam)

ax3 = plt.subplot(223, sharex=ax1)
paras = (numEval, numSite, numSam, 'PBC', 'true', flux)
ene = LoadEigVal(paras)
print(paras)
print(ene[0])
print(ene[1])
plt.plot(x[::inter], ene[::inter, 0]-ene[::inter, 0], '-o', color='red', label='0')
plt.plot(x[::inter], ene[::inter, 1]-ene[::inter, 0], '->', color='blue', label='$1$')
plt.plot(x[::inter], ene[::inter, 2]-ene[::inter, 0], '-d', color='green', label='$2$')
plt.plot(x[::inter], ene[::inter, 3]-ene[::inter, 0], '-x', color='yellow', label='$3$')
plt.plot(x[::inter], ene[::inter, 4]-ene[::inter, 0], '-h', color='black', label='$4$')
plt.xlabel('$J$')
plt.legend(loc='upper right', fontsize=8)
ax3.set_title('(c)')

numSite  =12
numSam = 500
x = np.linspace(step, step*numSam, num=numSam)

ax4 = plt.subplot(224, sharex=ax2)
paras = (numEval, numSite, numSam, 'PBC', 'true', flux)
ene = LoadEigVal(paras)
print(paras)
plt.plot(x[::inter], ene[::inter, 0]-ene[::inter, 0], '-o', color='red', label='0')
plt.plot(x[::inter], ene[::inter, 1]-ene[::inter, 0], '->', color='blue', label='$1$')
plt.plot(x[::inter], ene[::inter, 2]-ene[::inter, 0], '-d', color='green', label='$2$')
plt.plot(x[::inter], ene[::inter, 3]-ene[::inter, 0], '-x', color='yellow', label='$3$')
plt.plot(x[::inter], ene[::inter, 4]-ene[::inter, 0], '-h', color='black', label='$4$')
plt.xlabel('$J$')
plt.legend(loc='upper right', fontsize=8)
ax4.set_title('(d)')
  
#tit = '$t$-$J$ chain_varying_J_len_%s_step_%s_sam_%s_sigma_%s_con_%s.npy'
#plt.title(tit)
'''

image = 'energy_len_%s_step_%s_sam_%s_flux_%s_PI.pdf'
paras_image = (numSite, step, numSam, flux)
fig.tight_layout()
plt.savefig(image % paras_image, format='pdf')
plt.show()
