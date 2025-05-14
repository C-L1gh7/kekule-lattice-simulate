import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import AutoMinorLocator

t1 = 2.8
t2 = 0.5*t1
POE = 0.0
NOE = 0.0

main_L = 1.0  # length of main lattice
bond_L = 1.0  # bond length between lattices
modulus = 2 * main_L + bond_L
edgelength = 49.1 # edgelength = 4.1+3n(n=0,1,2,...)

kBT = 0.01
Mu = 0
Gamma = 0.001
h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength ** 2.
hS = -h_bar / S

mkpath="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)+""# 打印文件名

# set color
total = '#517556'
bb = '#F54C77'
cb = '#A09167'
eb = '#754F59'
ee = '#F4C64C'
ce = '#4C71F5'
cc = '#4CF563'

x = np.loadtxt(mkpath+'/h_omega.txt')
y_total = np.genfromtxt(mkpath + '/total_sigma_xx.txt', dtype=None)
y_B_B = np.genfromtxt(mkpath + '/B-B_sigma_xx.txt', dtype=None)
y_B_C = np.genfromtxt(mkpath + '/B-C_sigma_xx.txt', dtype=None)
y_B_E = np.genfromtxt(mkpath + '/B-E_sigma_xx.txt', dtype=None)
y_E_E = np.genfromtxt(mkpath + '/E-E_sigma_xx.txt', dtype=None)
y_E_C = np.genfromtxt(mkpath + '/E-C_sigma_xx.txt', dtype=None)
y_C_C = np.genfromtxt(mkpath + '/C-C_sigma_xx.txt', dtype=None)

fig, ax = plt.subplots(3, 1, sharex=True)
ax[0].plot(x, y_E_E, linestyle = '-', lw=0.5, color= ee, label = 'E-E')
ax[0].plot(x, y_B_C, linestyle = '--', lw=0.5, color= cb, label = 'C-B')
ax[0].plot(x, y_E_C, linestyle = ':', lw=0.5, color= ce, label = 'C-E')

ax[1].plot(x, y_C_C, linestyle = '-', lw=0.5, color= cc, label = 'C-C')
ax[1].plot(x, y_B_B, linestyle = '--', lw=0.5, color= bb, label = 'B-B')
ax[1].plot(x, y_B_E, linestyle = ':', lw=0.5, color= eb, label = 'E-B')

ax[2].plot(x, y_total, linestyle = '-', lw=0.5, color= total, label = 'Total')
ax[2].plot(x, y_B_B, linestyle = ':', lw=0.5, color= bb, label = 'B-B')
ax[2].plot(x, y_B_E, linestyle = '--', lw=0.5, color= eb, label = 'E-B')
ax[2].plot(x, y_B_C, linestyle = ':', lw=0.5, color= cb, label = 'C-B')
ax[2].plot(x, y_E_E, linestyle = '--', lw=0.5, color= ee, label = 'E-E')
ax[2].plot(x, y_E_C, linestyle = '-.', lw=0.5, color= ce, label = 'C-E')
ax[2].plot(x, y_C_C, linestyle = ':', lw=0.5, color= cc, label = 'C-C')

ax[0].set_xlim(0,2.5)
ax[1].set_xlim(0,2.5)


ax[0].legend(fontsize=7, markerscale=1, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
ax[1].legend(fontsize=7, markerscale=1, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
ax[2].legend(fontsize=7, markerscale=1, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left', ncol=2)

ax[0].minorticks_on() # 开启次刻度线
ax[1].minorticks_on()
ax[2].minorticks_on()

ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False,right=True, direction='in', width=0.5)
ax[1].tick_params(axis='both', which='both', top=True, labelbottom=False,right=True, direction='in', width=0.5)
ax[2].tick_params(axis='both', which='both', top=True, labelbottom=True,right=True, direction='in', width=0.5)

# 设置刻度线的长度和宽度
for a in ax:
    a.tick_params(axis='both', which='both', top=True, right=True, direction='in', width=0.5)
    a.tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
    a.tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
    a.xaxis.set_minor_locator(AutoMinorLocator(5))  # x轴主刻度之间有5个次刻度
    a.yaxis.set_minor_locator(AutoMinorLocator(5))  # y轴主刻度之间有5个次刻度

ax[0].text(-0.05, 0.95, '(a)', transform=ax[0].transAxes, fontsize=10, fontweight='bold')
ax[1].text(-0.05, 0.9, '(b)', transform=ax[1].transAxes, fontsize=10, fontweight='bold')
ax[2].text(-0.05, 0.9, '(c)', transform=ax[2].transAxes, fontsize=10, fontweight='bold')

ax[0].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')
ax[1].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')
ax[2].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')
ax[2].set_xlabel(r'$\hbar\omega\ [t]$')

plt.tight_layout()
plt.subplots_adjust(hspace=0)


plt.savefig(mkpath+'/optical conductivity 2.pdf')