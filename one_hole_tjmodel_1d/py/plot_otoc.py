#!/usr/bin/python3

import os
import math
import random
import numpy as np
from scipy import linalg as la

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid.inset_locator import inset_axes

beta = 20.0
timeNum = 500
timeStep = 0.02

# csFile = 'data/otoc_sz_J%s_beta%s_sigma%s_num%s_step%s.npy'
csFile = 'data/otoc_chargeSz36_J%s_beta%s_sigma%s_num%s_step%s.npy'

fs = 24 
start = 0
end = 120
inter = 2
fig = plt.figure(figsize=(8, 5))
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

x = np.linspace(0, timeNum*timeStep, num=timeNum)

ax1 = plt.subplot(111)
J = 0.3
npyParas = (J, beta, 'F', timeNum, timeStep)
cs0 = (np.load(csFile % npyParas))
npyParas = (J, beta, 'T', timeNum, timeStep)
cs1 = (np.load(csFile % npyParas))
# plt.plot(x[::inter], ee[::inter], '-o', color='red', label='$t$-$J$')
plt.plot(x[start:end:inter], cs0[start:end:inter], '.-', color='blue', label='$t$-$J$')
plt.plot(x[start:end:inter], cs1[start:end:inter], 'x-', color='red', label='$\sigma\cdot{t}$-$J$')
plt.xlabel('$T$', fontsize=fs)
# plt.ylabel('Commutator Square', fontsize=fs)
plt.legend(loc='best', frameon=False, fontsize=fs)
ax1.set_xlim([0, timeStep*(end-1)])
ax1.set_ylim(-0.002, 0.027)

xstart, xend = ax1.get_xlim()
# ax1.xaxis.set_ticks(np.arange(xstart, xend, 0.1))
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ystart, yend = ax1.get_ylim()
ax1.yaxis.set_ticks(np.arange(ystart+0.002, yend, 0.01))
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
ax1.tick_params(axis='both', labelsize=fs, direction='in')

image = 'otocSz_beta%s.pdf'
imageParas = (beta)
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
# plt.show()
