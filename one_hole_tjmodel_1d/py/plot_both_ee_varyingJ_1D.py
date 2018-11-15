import os
import math
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

fsize = 16

paras = 'False'
reader = csv.DictReader(open('./data/ee_varyingJ_size10_PBC_sigma%s.csv' % paras, 'r'), delimiter=',')

x = []
bEE0 = []
mEE0 = []
for row in reader:
    print(row)
    x.append(float(row['J']))
    bEE0.append(float(row['bee']))
    mEE0.append(float(row['mee']))

paras = 'True'
reader = csv.DictReader(open('./data/ee_varyingJ_size10_PBC_sigma%s.csv' % paras, 'r'), delimiter=',')

bEE1 = []
mEE1 = []
for row in reader:
    # print(row)
    bEE1.append(float(row['bee']))
    mEE1.append(float(row['mee']))

fig = plt.figure(figsize = (8, 4))
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
# plt.rc('font', family= 'serif')

inter = 10
ax1 = fig.add_subplot(121)
# ax1.set_title('(a) Etanglement Entropy', fontsize=fsize, x=0.3, y=0.9)
ax1.set_title('(a)', fontsize=fsize, x=0.1, y=0.9)
ax1.set_xlim(0.0, 40.0)
ax1.set_ylim(0.0, 3.5)
# Hide the right and top spines
# ax1.spines['right'].set_visible(False)
# ax1.spines['top'].set_visible(False)
ax1.axvline(x=17.1, color='green', linestyle='-.', linewidth=1)
ax1.axvline(x=28.5, color='green', linestyle='-.', linewidth=1)
ax1.plot(x[1::inter], bEE0[1::inter], '->', color='blue', label=r'$t$-$J$ Bipartite', linewidth=0.5)
ax1.plot(x[1::inter], bEE1[1::inter], '-x', color='red', label=r'$\sigma\cdot{t}$-$J$ Bipartite', linewidth=0.5)
ax1.plot(x[1::inter], mEE0[1::inter], '-o', color='green', label=r'$t$-$J$ Mutual', linewidth=0.5)
ax1.plot(x[1::inter], mEE1[1::inter], '-+', color='violet', label=r'$\sigma\cdot{t}$-$J$ Mutual', linewidth=0.5)
ax1.legend(loc='upper right', frameon=False, fontsize=fsize)
ax1.set_xlabel('$J/t$', fontsize=fsize)
# ax1.set_yticklabels('', visible=False)
# plt.setp(ax1.get_yticklabels(), visible=False)
# plt.ylabel('spin current aplitude', fontsize=fs)
# plt.setp(ax1.get_xticklabels(), fontsize=fs)
# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax1.tick_params(axis='both', labelsize=fsize, direction='in')

ystart, yend = ax1.get_ylim()
ax1.yaxis.set_ticks(np.arange(ystart, yend, 1.0))

paras = 'False'
reader = csv.DictReader(open('./data/eh_varyingJ_size10_PBC_sigma%s.csv' % paras, 'r'), delimiter=',')

x = []
bh0 = []
mh0 = []
for row in reader:
    x.append(float(row['J']))
    bh0.append(float(row['beh']))
    mh0.append(float(row['meh']))

paras = 'True'
reader = csv.DictReader(open('./data/eh_varyingJ_size10_PBC_sigma%s.csv' % paras, 'r'), delimiter=',')

bh1 = []
mh1 = []
for row in reader:
    bh1.append(float(row['beh']))
    mh1.append(float(row['meh']))

ax2 = fig.add_subplot(122)
ax2.set_title('(b)', fontsize=fsize, x=0.1, y=0.9)
ax2.set_xlim(0.0, 40.0)
ax2.set_ylim(0.0, 2.3)
ax2.axvline(x=17.1, color='green', linestyle='-.', linewidth=1)
ax2.axvline(x=28.5, color='green', linestyle='-.', linewidth=1)
ax2.plot(x[1::inter], bh0[1::inter], '->', color='blue', label=r'$t$-$J$ Bipartite', linewidth=0.5)
ax2.plot(x[1::inter], bh1[1::inter], '-x', color='red', label=r'$\sigma\cdot{t}$-$J$ Bipartite', linewidth=0.5)
ax2.plot(x[1::inter], mh0[1::inter], '-o', color='green', label=r'$t$-$J$ Mutual', linewidth=0.5)
ax2.plot(x[1::inter], mh1[1::inter], '-+', color='violet', label=r'$\sigma\cdot{t}$-$J$ Mutual', linewidth=0.5)
ax2.legend(loc='upper right', frameon=False, fontsize=fsize)
ax2.set_xlabel('$J/t$', fontsize=fsize)
# ax2.set_yticklabels('', visible=False)
# plt.setp(ax2.get_yticklabels(), visible=False)
ax2.tick_params(axis='both', labelsize=fsize, direction='in')

ystart, yend = ax2.get_ylim()
ax2.yaxis.set_ticks(np.arange(ystart, yend, 0.5))

image = 'both_ee_GS_size1D10_varyingJ.pdf'
# paras_image = (size, bounCon)
fig.tight_layout()
plt.savefig(image, format='PDF')
# plt.show()
