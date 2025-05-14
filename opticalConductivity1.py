import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import math
from tqdm import tqdm
from matplotlib.ticker import AutoMinorLocator

import basic_function as bf

edgelength1 = 48.1 # edgelength = 4.1+3n(n=0,1,2,...)
edgelength2 = 49.1 # edgelength = 4.1+3n(n=0,1,2,...)

kBT = 0.01
Mu = 0
Gamma = 0.001
h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength2 ** 2.
hS = -h_bar / S

####################################################################################

# color
with_corner_color = '#002450'
without_corner_color = '#E27E53'
Mu02 = '#B92636'
Mu047 = '#D5A351'

####################################################################################
mkpath1="result/anti-corner"+str(edgelength1)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
mkpath2="result/anti-corner"+str(edgelength2)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
mkpath3="result/anti-corner"+str(edgelength2)+"/mu=0.2"+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
mkpath4="result/anti-corner"+str(edgelength2)+"/mu=0.47"+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名

x1 = np.loadtxt(mkpath2+'/h_omega.txt')
y_nc = np.genfromtxt(mkpath1 + '/total_sigma_xx.txt', dtype=None) # with corner
y_c = np.genfromtxt(mkpath2 + '/total_sigma_xx.txt', dtype=None) # without corner
y_02 = np.genfromtxt(mkpath3 + '/total_sigma_xx.txt', dtype=None) # without corner
y_047 = np.genfromtxt(mkpath4 + '/total_sigma_xx.txt', dtype=None) # without corner

x2 = np.loadtxt(mkpath3+'/h_omega.txt')

fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)

ax[0].plot(x1, y_c, linestyle='-', lw=0.8, color=with_corner_color, label='L=49' )
ax[0].plot(x1, y_nc, linestyle = '--', lw=0.8, color=without_corner_color, label='L=48')
ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False,right=True, direction='in', width=0.5)
ax[0].set_ylim(bottom=-0.002)
ax[0].set_xlim(0, 2.5)
ax[0].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')
ax[0].legend(frameon=False)

ax[1].plot(x2, y_c, linestyle='-', lw=0.8, color=with_corner_color, label=r'$\mu = 0$', zorder=4)
ax[1].plot(x2, y_047, linestyle='-.', lw=0.8, color=Mu047, label=r'$\mu = 0.47$', zorder=4)
ax[1].plot(x2, y_02, linestyle=':', lw=0.8, color=Mu02, label=r'$\mu = 0.2$', zorder=4)
ax[1].set_xlabel(r'$\hbar\omega\ [t]$')
ax[1].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')
ax[1].tick_params(axis='both', which='both', top=True, labelbottom=True, right=True, direction='in', width=0.5)
ax[1].legend(frameon=False)
ax[1].set_ylim(bottom=-0.0025)
ax[1].set_xlim(0, 2.5)

ax[0].minorticks_on() # 开启次刻度线
ax[1].minorticks_on()
ax[0].tick_params(axis='both', which='major', length=3, width=0.4, direction='in') #设置主次刻度线的长度和宽度
ax[0].tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax[1].tick_params(axis='both', which='major', length=3, width=0.4, direction='in') #设置主次刻度线的长度和宽度
ax[1].tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))  # x轴主刻度之间有5个次刻度
ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))  # y轴主刻度之间有5个次刻度
ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))  # x轴主刻度之间有5个次刻度
ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))  # y轴主刻度之间有5个次刻度

ax[0].text(-0.1, 1.05, '(a)', transform=ax[0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')# 标注(a) (b)
ax[1].text(-0.1, 1.05, '(b)', transform=ax[1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

plt.subplots_adjust(hspace=0)
plt.savefig('E:/code/python/kekule_lattice/result/picture/optical conductivity 1.eps')
plt.clf()