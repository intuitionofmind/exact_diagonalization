import os
import math
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

fsize = 20

paras = 'False'
reader = csv.DictReader(open('./data/both_ee_GS_varyingJ_size52_bcOO_sigma%s.csv' % paras, 'r'), delimiter=',')

x = []
bEE0 = []
mEE0 = []
for row in reader:
    print(row)
    x.append(float(row['J']))
    bEE0.append(float(row['bee']))
    mEE0.append(float(row['mee']))

paras = 'True'
reader = csv.DictReader(open('./data/both_ee_GS_varyingJ_size52_bcOO_sigma%s.csv' % paras, 'r'), delimiter=',')

bEE1 = []
mEE1 = []
for row in reader:
    print(row)
    bEE1.append(float(row['bee']))
    mEE1.append(float(row['mee']))

fig = plt.figure(figsize = (8, 4))
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
# plt.rc('font', family= 'serif')

ax1 = fig.add_subplot(121)
ax1.set_title('(a) $t$-$J$ model', fontsize=fsize, x=0.2, y=0.9)
ax1.set_xlim(0.0, 10.0)
ax1.set_ylim(0.5, 5.0)
# ax.axvline(x=2.05, color='red', linestyle='-.', linewidth=1)
ax1.plot(x[1::2], bEE0[1::2], '->', color='blue', label='bEE', linewidth=0.5)
ax1.plot(x[1::2], mEE0[1::2], '-o', color='red', label='mEE', linewidth=0.5)
ax1.legend(loc='upper right', frameon=False, fontsize=fsize)
ax1.set_xlabel('$J/t$', fontsize=fsize)
# ax1.set_yticklabels('', visible=False)
# plt.setp(ax1.get_yticklabels(), visible=False)
# plt.ylabel('spin current aplitude', fontsize=fs)
# plt.setp(ax1.get_xticklabels(), fontsize=fs)
# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax1.tick_params(axis='both', labelsize=fsize, direction='in')
# ax1.tick_params('x', labelsize=16)
# ax1.tick_params('y', labelsize=16)

ax2 = fig.add_subplot(122, sharey=ax1)
ax2.set_title('(b) $\sigma\cdot{t}$-$J$ model', fontsize=fsize, x=0.25, y=0.9)
ax2.set_xlim(0.0, 10.0)
ax2.set_ylim(0.0, 5.0)
ax2.plot(x[1::2], bEE1[1::2], '->', color='blue', label='bEE', linewidth=0.5)
ax2.plot(x[1::2], mEE1[1::2], '-o', color='red', label='mEE', linewidth=0.5)
ax2.legend(loc='upper right', frameon=False, fontsize=fsize)
ax2.set_xlabel('$J/t$', fontsize=fsize)
# ax2.set_yticklabels('', visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.tick_params(axis='both', labelsize=fsize, direction='in')

image = 'both_ee_GS_size52_varyingJ.pdf'
# paras_image = (size, bounCon)
fig.tight_layout()
plt.savefig(image, format='PDF')
# plt.show()
