#!/usr/bin/python3

import os
import csv
import numpy as np
from scipy import linalg as la
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from auxi_1d import numSite, numEval, cut, fs

level = 299
numTime = 1000
inter = 4
numPlot = int(numTime/inter)

plotInter = 2
fs = 24

fig = plt.figure(figsize=(10, 10))
mpl.rcParams['axes.linewidth'] = 2.0
plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

deltaT = 0.1
ax1 = plt.subplot(211)

csvFile = 'data/test_time_evolution_size52_level0_J5.3-0.3.csv'
reader = csv.DictReader(open(csvFile, 'r'), delimiter=',')
t = []
ee = []
ee0 = []
ee1 = []
for row in reader:
    t.append(float(row['t']))
    ee.append(float(row['eeFull']))
    ee0.append(float(row['eePartial0']))
    ee1.append(float(row['eePartial1']))

ax1.set_title('(a) $t$-$J$, $l=0, J/t=5.3\longrightarrow{0.3}$', fontsize=fs, x=0.5, y=0.1)
plt.plot(t[::plotInter], ee[::plotInter], '.-', color='blue', label='bEE-Full')
plt.plot(t[::plotInter], ee0[::plotInter], 'x-', color='red', label='bEE-Partial-200')
plt.plot(t[::plotInter], ee1[::plotInter], 'o-', color='green', label='bEE-Partial-630')
plt.legend(loc='best', frameon=False, fontsize=fs)
ax1.set_xlim([0, 10.0])
ax1.set_ylim([2.6, 3.8])
# plt.yticks(yt, fontsize=fs)
xstart, xend = ax1.get_xlim()
ax1.xaxis.set_ticks(np.arange(xstart, xend, 2.0))
ystart, yend = ax1.get_ylim()
ax1.yaxis.set_ticks(np.arange(ystart, yend, 0.3))
ax1.tick_params(axis='both', labelsize=fs, direction='in')


ax2 = plt.subplot(212)

csvFile = 'data/test_time_evolution_size52_level100_J5.3-0.3.csv'
reader = csv.DictReader(open(csvFile, 'r'), delimiter=',')
t = []
ee = []
ee0 = []
ee1 = []
for row in reader:
    t.append(float(row['t']))
    ee.append(float(row['eeFull']))
    ee0.append(float(row['eePartial0']))
    ee1.append(float(row['eePartial1']))

ax2.set_title('(b) $t$-$J$, $l=100, J/t=5.3\longrightarrow{0.3}$', fontsize=fs, x=0.5, y=0.1)
plt.plot(t[::plotInter], ee[::plotInter], '.-', color='blue', label='bEE-Full')
plt.plot(t[::plotInter], ee0[::plotInter], 'x-', color='red', label='bEE-Partial-200')
plt.plot(t[::plotInter], ee1[::plotInter], 'o-', color='green', label='bEE-Partial-630')
plt.legend(loc='best', frameon=False, fontsize=fs)
ax2.set_xlim([0, 10.0])
ax2.set_ylim([3.2, 3.8])
xstart, xend = ax2.get_xlim()
ax2.xaxis.set_ticks(np.arange(xstart, xend, 2.0))
ystart, yend = ax2.get_ylim()
ax2.yaxis.set_ticks(np.arange(ystart, yend, 0.3))
ax2.tick_params(axis='both', labelsize=fs, direction='in')

image = 'test_time_evolution.pdf'
# imageParas = (numSite, 'P')
fig.tight_layout()
plt.savefig(image, format='pdf')
# plt.show()
