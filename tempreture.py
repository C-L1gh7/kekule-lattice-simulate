# Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from tqdm import tqdm

import basic_function as bf

edgelength = 49.1

kBT_1 = 0.01
Mu_1 = 0
Gamma_1 = 0.0005

kBT_2 = 0.025
Mu_2 = 0
Gamma_2 = 0.0005

kBT_3 = 0.02
Mu_3 = 0
Gamma_3 = 0.0005

kBT_4 = 0.003
Mu_4 = 0
Gamma_4 = 0.0005

kBT_5 = 0.001
Mu_5 = 0
Gamma_5 = 0.0005

mkpath_1="result/anti-corner"+str(edgelength)+"/mu="+str(Mu_1)+"kBT="+str(kBT_1)+"gamma="+str(Gamma_1)# 打印文件名
mkpath_2="result/anti-corner"+str(edgelength)+"/mu="+str(Mu_2)+"kBT="+str(kBT_2)+"gamma="+str(Gamma_2)# 打印文件名
mkpath_3="result/anti-corner"+str(edgelength)+"/mu="+str(Mu_3)+"kBT="+str(kBT_3)+"gamma="+str(Gamma_3)# 打印文件名
mkpath_4="result/anti-corner"+str(edgelength)+"/mu="+str(Mu_4)+"kBT="+str(kBT_4)+"gamma="+str(Gamma_4)# 打印文件名
mkpath_5="result/anti-corner"+str(edgelength)+"/mu="+str(Mu_5)+"kBT="+str(kBT_5)+"gamma="+str(Gamma_5)# 打印文件名

x = np.loadtxt(mkpath_1+'/h_omega.txt')
y1 = np.genfromtxt(mkpath_1 + '/total_sigma_xx.txt', dtype=None)
y2 = np.genfromtxt(mkpath_2 + '/total_sigma_xx.txt', dtype=None)
y3 = np.genfromtxt(mkpath_3 + '/total_sigma_xx.txt', dtype=None)
y4 = np.genfromtxt(mkpath_4 + '/total_sigma_xx.txt', dtype=None)
y5 = np.genfromtxt(mkpath_5 + '/total_sigma_xx.txt', dtype=None)

plt.plot(x, y1, color= '#F553DA', linestyle='-', lw=0.8, label = r'$T = 115.9K$')
plt.plot(x, y2, color= '#53A8F5', linestyle='-.', lw=0.8, label = r'$T = 289.7K$')
plt.plot(x, y3, color= '#756372',linestyle=':', lw=0.8, label = r'$T = 231.8K$')
plt.plot(x, y4, color = '#F5B753', linestyle='--', lw=0.8, label = r'$T = 34.8K$')
plt.plot(x, y5, color = '#79F553', linestyle='-.', lw=0.8, label = r'$T = 11.6K$')

plt.legend(frameon=False) # 显示图例
plt.xlabel(r'$\hbar\omega\ [t]$', fontsize=14)
plt.ylabel(r'$\sigma_{xx}\ [{\tilde t}^2 e^2/\hbar]$', fontsize=14)
plt.savefig('result/T_sigma_xx_total.eps')  # 保存为eps文件
plt.clf

y1_BB = np.genfromtxt(mkpath_1 + '/B-B_sigma_xx.txt', dtype=None)
y1_BC = np.genfromtxt(mkpath_1 + '/B-C_sigma_xx.txt', dtype=None)
y1_CC = np.genfromtxt(mkpath_1 + '/C-C_sigma_xx.txt', dtype=None)

y2_BB = np.genfromtxt(mkpath_2 + '/B-B_sigma_xx.txt', dtype=None)
y2_BC = np.genfromtxt(mkpath_2 + '/B-C_sigma_xx.txt', dtype=None)
y2_CC = np.genfromtxt(mkpath_2 + '/C-C_sigma_xx.txt', dtype=None)

y3_BB = np.genfromtxt(mkpath_3 + '/B-B_sigma_xx.txt', dtype=None)
y3_BC = np.genfromtxt(mkpath_3 + '/B-C_sigma_xx.txt', dtype=None)
y3_CC = np.genfromtxt(mkpath_3 + '/C-C_sigma_xx.txt', dtype=None)

y4_BB = np.genfromtxt(mkpath_4 + '/B-B_sigma_xx.txt', dtype=None)
y4_BC = np.genfromtxt(mkpath_4 + '/B-C_sigma_xx.txt', dtype=None)
y4_CC = np.genfromtxt(mkpath_4 + '/C-C_sigma_xx.txt', dtype=None)

y5_BB = np.genfromtxt(mkpath_5 + '/B-B_sigma_xx.txt', dtype=None)
y5_BC = np.genfromtxt(mkpath_5 + '/B-C_sigma_xx.txt', dtype=None)
y5_CC = np.genfromtxt(mkpath_5 + '/C-C_sigma_xx.txt', dtype=None)

data = [
    {'total':y1,'BB':y1_BB,'BC':y1_BC,'CC':y1_CC},
    {'total':y2,'BB':y2_BB,'BC':y2_BC,'CC':y2_CC},
    {'total':y3,'BB':y3_BB,'BC':y3_BC,'CC':y3_CC},
    {'total':y4,'BB':y4_BB,'BC':y4_BC,'CC':y4_CC},
    {'total':y5,'BB':y5_BB,'BC':y5_BC,'CC':y5_CC}
]

y_max = max(max(data[i][key]) for i in range(len(data)) for key in data[i])
fig, axs = plt.subplots(5, 1,sharex=True, sharey=True)# 合并x，y轴
axs = axs.ravel()
for i, ax in enumerate(axs):
    ax.plot(x, data[i]['total'], label='Total', color=(55/255,103/255,149/255), linestyle='-', lw=0.8)
    ax.plot(x, data[i]['BB'], label='B-B', color=(114/255,188/255,213/255), linestyle='--', lw=0.8)
    ax.plot(x, data[i]['CC'], label='C-C', color=(255/255,208/255,111/255), linestyle='-.', lw=0.8)
    ax.plot(x, data[i]['BC'], label='B-C', color=(231/255,98/255,84/255), linestyle=':', lw=0.8)
    ax.set_aspect(aspect=50.0)
    # ax.set_xlim(0, 2)  # 设置x轴范围
    # ax.set_ylim(0, y_max)    # 设置y轴范围
handles, labels = axs[0].get_legend_handles_labels()
plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05, hspace=0.0)
fig.legend(handles, labels,frameon=False) # 显示图例)
plt.savefig('result/T_sigma_xx_all.eps')  # 保存为eps文件