import os
import math
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

fsize = 20
inter = 1

paras = 'False'
# reader = csv.DictReader(open('./data/both_ee_excited_size44_bcPP_sigma%s.csv' % paras, 'r'), delimiter=',')
reader = csv.DictReader(open('./data/both_ee_excited_J0.3_size52_bcOO_sigma%s.csv' % paras, 'r'), delimiter=',')

x0 = []
bEE0 = []
mEE0 = []
for row in reader:
    x0.append(float(row['ene_density']))
    bEE0.append(float(row['bee']))
    mEE0.append(float(row['mee']))

paras = 'True'
# reader = csv.DictReader(open('./data/both_ee_excited_size44_bcPP_sigma%s.csv' % paras, 'r'), delimiter=',')
reader = csv.DictReader(open('./data/both_ee_excited_J0.3_size52_bcOO_sigma%s.csv' % paras, 'r'), delimiter=',')

x1 = []
bEE1 = []
mEE1 = []
for row in reader:
    x1.append(float(row['ene_density']))
    bEE1.append(float(row['bee']))
    mEE1.append(float(row['mee']))

paras = 'False'
# reader = csv.DictReader(open('./data/both_ee_excited_size44_bcPP_sigma%s.csv' % paras, 'r'), delimiter=',')
reader = csv.DictReader(open('./data/both_ee_excited_J40.3_size52_bcOO_sigma%s.csv' % paras, 'r'), delimiter=',')

x2 = []
bEE2 = []
mEE2 = []
for row in reader:
    x2.append(float(row['ene_density']))
    bEE2.append(float(row['bee']))
    mEE2.append(float(row['mee']))

paras = 'True'
# reader = csv.DictReader(open('./data/both_ee_excited_size44_bcPP_sigma%s.csv' % paras, 'r'), delimiter=',')
reader = csv.DictReader(open('./data/both_ee_excited_J40.3_size52_bcOO_sigma%s.csv' % paras, 'r'), delimiter=',')

x3 = []
bEE3 = []
mEE3 = []
for row in reader:
    x3.append(float(row['ene_density']))
    bEE3.append(float(row['bee']))
    mEE3.append(float(row['mee']))

print('Loading finished.')

fig = plt.figure(figsize = (8, 8))
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
# plt.rc('font', family= 'serif')

ax1 = fig.add_subplot(221)
ax1.set_title('(a) $t$-$J$, J/t=0.3', fontsize=fsize, x=0.3, y=0.9)
ax1.set_xlim(x0[0], x0[999])
ax1.set_ylim(0.0, 5.0)
# ax.axvline(x=2.05, color='red', linestyle='-.', linewidth=1)
ax1.plot(x0[::inter], bEE0[::inter], '.', color='blue', label='bEE', linewidth=0.5)
ax1.plot(x0[::inter], mEE0[::inter], 'x', color='red', label='mEE', linewidth=0.5)
ax1.legend(loc='upper right', frameon=False, fontsize=fsize)
# ax1.set_xlabel('$E/L$', fontsize=fsize)
ax1.tick_params(axis='both', labelsize=fsize, direction='in')

xstart, xend = ax1.get_xlim()
ax1.xaxis.set_ticks(np.arange(xstart, xend, 0.15))
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
ystart, yend = ax1.get_ylim()
ax1.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))


ax2 = fig.add_subplot(222, sharey=ax1)
ax2.set_title('(b) $\sigma\cdot{t}$-$J$, J/t=0.3', fontsize=fsize, x=0.35, y=0.9)
ax2.set_xlim(x1[0], x1[999])
ax2.set_ylim(0.0, 5.0)
ax2.plot(x1[::inter], bEE1[::inter], '.', color='blue', label='bEE', linewidth=0.5)
ax2.plot(x1[::inter], mEE1[::inter], 'x', color='red', label='mEE', linewidth=0.5)
ax2.legend(loc='upper right', frameon=False, fontsize=fsize)
# ax2.set_xlabel('$E/L$', fontsize=fsize)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.tick_params(axis='both', labelsize=fsize, direction='in')

xstart, xend = ax2.get_xlim()
ax2.xaxis.set_ticks(np.arange(xstart, xend, 0.15))
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
ystart, yend = ax2.get_ylim()
ax2.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))


ax3 = fig.add_subplot(223)
ax3.set_title('(c) $t$-$J$, J/t=40.3', fontsize=fsize, x=0.3, y=0.9)
ax3.set_xlim(x2[0], x2[999])
ax3.set_ylim(0.0, 5.0)
# ax.axvline(x=2.05, color='red', linestyle='-.', linewidth=1)
ax3.plot(x2[::inter], bEE2[::inter], '.', color='blue', label='bEE', linewidth=0.5)
ax3.plot(x2[::inter], mEE2[::inter], 'x', color='red', label='mEE', linewidth=0.5)
ax3.legend(loc='upper right', frameon=False, fontsize=fsize)
ax3.set_xlabel('$E/L$', fontsize=fsize)
ax3.tick_params(axis='both', labelsize=fsize, direction='in')

xstart, xend = ax3.get_xlim()
ax3.xaxis.set_ticks(np.arange(xstart, xend, 6.0))
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
ystart, yend = ax3.get_ylim()
ax3.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

ax4 = fig.add_subplot(224, sharey=ax3)
ax4.set_title('(d) $\sigma\cdot{t}$-$J$, J/t=40.3', fontsize=fsize, x=0.35, y=0.9)
ax4.set_xlim(x3[0], x3[999])
ax4.set_ylim(0.0, 5.0)
ax4.plot(x3[::inter], bEE3[::inter], '.', color='blue', label='bEE', linewidth=0.5)
ax4.plot(x3[::inter], mEE3[::inter], 'x', color='red', label='mEE', linewidth=0.5)
ax4.legend(loc='upper right', frameon=False, fontsize=fsize)
ax4.set_xlabel('$E/L$', fontsize=fsize)
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.tick_params(axis='both', labelsize=fsize, direction='in')

xstart, xend = ax4.get_xlim()
ax4.xaxis.set_ticks(np.arange(xstart, xend, 6.0))
ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
ystart, yend = ax4.get_ylim()
ax4.yaxis.set_ticks(np.arange(ystart+1.0, yend, 1.0))
ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))


image = 'both_ee_excited_size52_bcOO.pdf'
# paras_image = (size, bounCon)
fig.tight_layout()
plt.savefig(image, format='PDF')
# plt.show()
