# Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from tqdm import tqdm

import basic_function as bf

####################################################################################

edgelength = 48.1
kBT = 0.01
Mu = 0
Gamma = 0.001

####################################################################################

#设置存储位置
mkpath="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)+"/yy"# 打印文件名

x = np.loadtxt(mkpath+'/h_omega.txt')
y_total = np.genfromtxt(mkpath + '/total_sigma_xx.txt', dtype=None)


plt.plot(x, y_total, color= (55/255,103/255,149/255), linestyle='-', lw=0.8, label = 'Total')
# 添加标题和轴标签
plt.legend(frameon=False) # 显示图例
plt.xlabel(r'$\hbar\omega\ [t]$', fontsize=14)
plt.ylabel(r'$\sigma_{xx}\ [{\tilde t}^2 e^2/\hbar]$', fontsize=14)
plt.savefig(mkpath+'/sigma_yy.pdf')  # 保存为eps文件
