# Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from tqdm import tqdm
from numba import njit
from scipy.sparse import csr_matrix

import basic_function as bf

data1 = {
    'edgelength': 49.1,
    'kBT': 0.01,
    'Mu': 0,
    'Gamma': 0.001}

data2 = {
    'edgelength': 49.1,
    'kBT': 0.01,
    'Mu': 0,
    'Gamma': 0.001}

mkpath1="result/anti-corner"+str(data1['edgelength'])+"/mu="+str(data1['Mu'])+"kBT="+str(data1['kBT'])+"gamma="+str(data1['Gamma'])# 打印文件名
mkpath2="result/anti-corner"+str(data2['edgelength'])+"/mu="+str(data2['Mu'])+"kBT="+str(data2['kBT'])+"gamma="+str(data2['Gamma'])+'/yy'# 打印文件名
x1 = np.loadtxt(mkpath1+'/h_omega.txt')
y1 = np.genfromtxt(mkpath1 + '/total_sigma_xx.txt', dtype=None)

x2 = np.loadtxt(mkpath1+'/h_omega.txt')
y2 = np.genfromtxt(mkpath2 + '/total_sigma_xx.txt', dtype=None)

plt.plot(x1, y1, color= 'red', linestyle='-', lw=0.5, label = 'corner')
plt.plot(x2, y2, color= 'blue', linestyle='--', lw=0.5, label = 'no corner')

plt.legend(frameon=False) # 显示图例
plt.xlabel(r'$\hbar\omega\ [t]$', fontsize=14)
plt.ylabel(r'$\sigma_{xx}\ [{\tilde t}^2 e^2/\hbar]$', fontsize=14)
plt.savefig('sigma_test.pdf')  # 保存为eps文件