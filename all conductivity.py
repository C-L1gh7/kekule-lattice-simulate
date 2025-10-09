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
edgelength1 = 48.1 # edgelength = 4.1+3n(n=0,1,2,...)
edgelength2 = 49.1 # edgelength = 4.1+3n(n=0,1,2,...)

kBT = 0.01
Mu = 0
Gamma = 0.001
h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength ** 2.
hS = -h_bar / S

mkpath="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)+""# 打印文件名
mkpath1="result/anti-corner"+str(edgelength1)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
mkpath2="result/anti-corner"+str(edgelength2)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名

# set color
with_corner_color_Between = "#FFB7F9"
without_corner_color_Between = "#FFA396"

with_corner_color = "#FF6CF3"
without_corner_color = "#FF6652"

total = '#517556'
bb_Between = "#8AC7FF"
cb_Between = "#F2DBFD"
eb_Between = "#92FFF0"
ee_Between = "#9DFFBA"
ce_Between = "#FFFD72"
cc_Between = "#FF9C9C"

bb = "#40A3FF"
cb = "#CC5FFF"
eb = "#24D8C0"
ee = "#1BA846"
ce = "#D6C800"
cc = "#EC5252"

alpha=0.3
lw=0.3

x = np.loadtxt(mkpath+'/h_omega.txt')
y_total = np.genfromtxt(mkpath + '/total_sigma_xx.txt', dtype=None)
y_B_B = np.genfromtxt(mkpath + '/B-B_sigma_xx.txt', dtype=None)
y_B_C = np.genfromtxt(mkpath + '/B-C_sigma_xx.txt', dtype=None)
y_B_E = np.genfromtxt(mkpath + '/B-E_sigma_xx.txt', dtype=None)
y_E_E = np.genfromtxt(mkpath + '/E-E_sigma_xx.txt', dtype=None)
y_E_C = np.genfromtxt(mkpath + '/E-C_sigma_xx.txt', dtype=None)
y_C_C = np.genfromtxt(mkpath + '/C-C_sigma_xx.txt', dtype=None)
y_nc = np.genfromtxt(mkpath1 + '/total_sigma_xx.txt', dtype=None) # without corner
y_c = np.genfromtxt(mkpath2 + '/total_sigma_xx.txt', dtype=None) # with corner

fig, ax = plt.subplots(3, 1, sharex=True)

# 统一的legend字体大小
legend_fontsize = 8

ax[0].plot(x, y_c, linestyle='--', lw=lw, color=with_corner_color, label='L=49' )
ax[0].fill_between(x, y_c, color=with_corner_color_Between, alpha=alpha)
ax[0].plot(x, y_nc, linestyle = ':', lw=lw, color=without_corner_color, label='L=48')
ax[0].fill_between(x, y_nc, color=without_corner_color_Between, alpha=alpha)
ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False,right=True, direction='in', width=0.5)
# 修复：移除手动设置的ylim，让matplotlib自动调整
# ax[0].set_ylim(bottom=-0.002)
ax[0].set_xlim(0, 2.2)
ax[0].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')
ax[0].legend(frameon=False, loc='upper left', fontsize=legend_fontsize)

ax[1].plot(x, y_B_B, linestyle = '--', lw=lw, color= bb, label = 'B-B')
ax[1].fill_between(x, y_B_B, color=bb_Between, alpha=alpha)
ax[1].plot(x, y_B_E, linestyle = ':', lw=lw, color= eb, label = 'E-B')
ax[1].fill_between(x, y_B_E, color=eb_Between, alpha=alpha)
# ax[1].plot(x, y_C_C, linestyle = '-', lw=lw, color= cc, label = 'C-C')
# ax[1].fill_between(x, y_C_C, color=cc_Between, alpha=alpha)
ax[1].minorticks_on() # 开启次刻度线
ax[1].set_xlim(0,2.2)
ax[1].legend(fontsize=legend_fontsize, markerscale=1, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
ax[1].tick_params(axis='both', which='both', top=True, labelbottom=False,right=True, direction='in', width=0.5)

ax[2].plot(x, y_B_C, linestyle = ':', lw=lw, color= cb, label = 'C-B')
ax[2].fill_between(x, y_B_C, color=cb_Between, alpha=alpha)
ax[2].plot(x, y_E_C, linestyle = ':', lw=lw, color= ce, label = 'C-E')
ax[2].fill_between(x, y_E_C, color=ce_Between, alpha=alpha)
ax[2].plot(x, y_E_E, linestyle = '-', lw=lw, color= ee, label = 'E-E')
ax[2].fill_between(x, y_E_E, color=ee_Between, alpha=alpha)
ax[2].minorticks_on() # 开启次刻度线
ax[2].set_xlim(0,2.2)
ax[2].set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.2])
ax[2].legend(fontsize=legend_fontsize, markerscale=1, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
ax[2].tick_params(axis='both', which='both', top=True, labelbottom=True,right=True, direction='in', width=0.5)
ax[2].set_yticks([0.000, 0.002])
ax[2].set_yticklabels(['0.000', '0.002'])

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
ax[2].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$', fontsize=14)
ax[2].set_xlabel(r'$\hbar\omega\ [t]$')

plt.tight_layout()
plt.subplots_adjust(hspace=0)

plt.savefig('result/picture/optical conductivity 2.svg')
# plt.show()