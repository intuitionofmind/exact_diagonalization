#!/usr/bin/python3

import os
import csv
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from auxi_1d import numSite, numEval, cut, fs

level = 299
numTime = 1000
inter = 4
numPlot = int(numTime/inter)

plotInter = 5

dataDir = '/Users/wayne/Downloads/data/'
csvFile = 'timeEE_J%sto%s_step%s_state%s_bc%s_sigma%s.csv'

fig = plt.figure(figsize=(6, 6))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

deltaT = 0.1
x = np.linspace(deltaT, deltaT*numTime, num=numTime-1)
xl = np.log(x)
ax1 = plt.subplot(211)
csvParas = (40.3, 0.3, deltaT, level, 'P', 'T')
reader = csv.DictReader(open(os.path.join(dataDir, csvFile) % csvParas, 'r'), delimiter=',')
mee = []
bee = []
for row in reader:
    mee.append(float(row['mee']))
    bee.append(float(row['see']))
ax1.set_title('(a) $t$-$J$, $J=0.3$', fontsize=fs, x=0.5, y=0.1)
plt.plot(xl[::plotInter], mee[1::plotInter], '.-', color='red', label='mEE')
plt.plot(xl[::plotInter], bee[1::plotInter], 'x-', color='blue', label='bEE')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
# plt.setp(ax1.get_xticklabels(), visible=False, fontsize=fs)
plt.xticks(np.log([50, 100]), [50, 100], fontsize=fs)
ax1.set_xlim([np.log(deltaT), np.log(deltaT*numTime)])
ax1.set_ylim([0, 5])
# plt.yticks(yt, fontsize=fs)
print(csvParas)
print(np.mean(mee[int(numTime/2):]), np.std(mee[int(numTime/2):])/(np.sqrt(numTime/2)))
print(np.mean(bee[int(numTime/2):]), np.std(bee[int(numTime/2):])/(np.sqrt(numTime/2)))

deltaT = 0.0025
x = np.linspace(deltaT, deltaT*numTime, num=numTime-1)
xl = np.log(x)
ax2 = plt.subplot(212)
csvParas = (0.3, 40.3, deltaT, level, 'P', 'T')
reader = csv.DictReader(open(os.path.join(dataDir, csvFile) % csvParas, 'r'), delimiter=',')
mee = []
bee = []
for row in reader:
    mee.append(float(row['mee']))
    bee.append(float(row['see']))
ax2.set_title('(b) $t$-$J$, $J=40.3$', fontsize=fs, x=0.5, y=0.1)
plt.plot(xl[::plotInter], mee[1::plotInter], '.-', color='red', label='mEE')
plt.plot(xl[::plotInter], bee[1::plotInter], 'x-', color='blue', label='bEE')
plt.xticks(np.log([1.25, 2.5]), [1.25, 2.5], fontsize=fs)
ax2.set_xlim([np.log(deltaT), np.log(deltaT*numTime)])
ax2.set_ylim([0, 5])
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.xlabel('$Tt$', fontsize=fs)
print(csvParas)
print(np.mean(mee[int(numTime/2):]), np.std(mee[int(numTime/2):])/(np.sqrt(numTime/2)))
print(np.mean(bee[int(numTime/2):]), np.std(bee[int(numTime/2):])/(np.sqrt(numTime/2)))

'''
ax4 = plt.subplot(224, sharey=ax3)
csvParas = (0.3, 40.3, deltaT, level, 'P', 'T')
reader = csv.DictReader(open(os.path.join(dataDir, csvFile) % csvParas, 'r'), delimiter=',')
mee = []
bee = []
for row in reader:
    mee.append(float(row['mee']))
    bee.append(float(row['see']))
ax4.set_title('(d) $\sigma\cdot{t}$-$J$, $J=40.3$', fontsize=fs, x=0.25, y=0.1)
plt.plot(x[::plotInter], mee[::plotInter], '.-', color='red', label='mEE')
plt.plot(x[::plotInter], bee[::plotInter], 'x-', color='blue', label='bEE')
plt.legend(loc='lower right', frameon=False, fontsize=fs)
plt.setp(ax4.get_yticklabels(), visible=False, fontsize=fs)
plt.xlabel('$Tt$', fontsize=fs)
# plt.xticks(xt, fontsize=fs)
ax4.set_xlim([0, deltaT*numTime/2])
ax4.set_ylim([0, 5])
print(csvParas)
print(np.mean(mee[int(numTime/2):]), np.std(mee[int(numTime/2):])/(np.sqrt(numTime/2)))
print(np.mean(bee[int(numTime/2):]), np.std(bee[int(numTime/2):])/(np.sqrt(numTime/2)))
'''

image = 'timeEE_len%s_bc%s.pdf'
imageParas = (numSite, 'P')
fig.tight_layout()
plt.savefig(image % imageParas, format='pdf')
# plt.show()
