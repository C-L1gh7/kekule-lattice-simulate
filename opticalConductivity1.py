import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import math
from tqdm import tqdm
from matplotlib.ticker import AutoMinorLocator

import basic_function as bf
####################################################################################
edgelength = 49.1 # edgelength = 4.1+3n(n=0,1,2,...)

kBT = 0.01
Mu = 0
Gamma = 0.001
h_bar = 1.0

####################################################################################

alpha=0.3
lw=0.8
legend_lw = 1.2  # legend线条粗细
legend_fontsize = 14  # legend字体大小

mkpath0="result/anti-corner"+str(edgelength)+"/mu=0"+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
mkpath02="result/anti-corner"+str(edgelength)+"/mu=0.2"+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
mkpath047="result/anti-corner"+str(edgelength)+"/mu=0.47"+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名

mkpath_1="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(0.0005)# 打印文件名
mkpath_2="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(0.001)# 打印文件名
mkpath_3="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(0.005)# 打印文件名
# mkpath_4="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(0.01)# 打印文件名
# mkpath_5="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(0.02)# 打印文件名
mkpath_6="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(0.03)# 打印文件名

x_mu = np.loadtxt(mkpath0+'/h_omega.txt')
y_0 = np.genfromtxt(mkpath0 + '/total_sigma_xx.txt', dtype=None)
y_02 = np.genfromtxt(mkpath02 + '/total_sigma_xx.txt', dtype=None)
y_047 = np.genfromtxt(mkpath047 + '/total_sigma_xx.txt', dtype=None)

x1 = np.loadtxt(mkpath_1+'/h_omega.txt')
y1 = np.genfromtxt(mkpath_1 + '/total_sigma_xx.txt', dtype=None)
x2 = np.loadtxt(mkpath_2+'/h_omega.txt')
y2 = np.genfromtxt(mkpath_2 + '/total_sigma_xx.txt', dtype=None)
x3 = np.loadtxt(mkpath_3+'/h_omega.txt')
y3 = np.genfromtxt(mkpath_3 + '/total_sigma_xx.txt', dtype=None)
# x4 = np.loadtxt(mkpath_4+'/h_omega.txt')
# y4 = np.genfromtxt(mkpath_4 + '/total_sigma_xx.txt', dtype=None)
# x5 = np.loadtxt(mkpath_5+'/h_omega.txt')
# y5 = np.genfromtxt(mkpath_5 + '/total_sigma_xx.txt', dtype=None)
x6 = np.loadtxt(mkpath_6+'/h_omega.txt')
y6 = np.genfromtxt(mkpath_6 + '/total_sigma_xx.txt', dtype=None)

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 8))

# 图a绘制
ax[0].plot(x_mu, y_047, linestyle='-', lw=lw, color="#7EB833", label=r'$\mu = 0.47t_1$', zorder=4)
# ax[0].fill_between(x_mu, y_047, color="#87E2E6", alpha=alpha)
ax[0].plot(x_mu, y_02, linestyle='-', lw=lw, color="#DD58D6", label=r'$\mu = 0.2t_1$', zorder=4)
# ax[0].fill_between(x_mu, y_02, color="#FF8FDD", alpha=alpha)
ax[0].plot(x_mu, y_0, linestyle='--', lw=lw, color="#18A7AC", label=r'$\mu = 0t_1$', zorder=4)
# ax[0].fill_between(x_mu, y_0, color="#B0FF91", alpha=alpha)

# 计算图a实际绘制数据的y轴最大值（在x范围[0,2]内）
x_mask_a = (x_mu >= 0) & (x_mu <= 2)
y_plotted_a = np.concatenate([y_047[x_mask_a], y_02[x_mask_a], y_0[x_mask_a]])
y_min_a = np.min(y_plotted_a)
y_max_a = np.max(y_plotted_a)

ax[0].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$', fontsize=16)
ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False, right=True, direction='in', width=0.5)  # 关闭x轴刻度标签

# 设置legend到左上，加粗图例线条
# 获取当前图例的handles和labels
handles, labels = ax[0].get_legend_handles_labels()
# 反转顺序
handles = handles[::-1]
labels = labels[::-1]
# 重新生成legend（顺序已反转）
leg0 = ax[0].legend(handles, labels, handlelength=2, frameon=False,
                    loc='upper left', fontsize=legend_fontsize)
# 设置图例线条宽度
for legline in leg0.get_lines():
    legline.set_linewidth(legend_lw)

ax[0].set_ylim(-0.0005, y_max_a * 1.1)
ax[0].set_xlim(0, 2.2)
ax[0].minorticks_on() # 开启次刻度线
ax[0].tick_params(axis='both', which='major', length=3, width=0.4, direction='in') #设置主次刻度线的长度和宽度
ax[0].tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))  # x轴主刻度之间有5个次刻度
ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))  # y轴主刻度之间有5个次刻度
ax[0].text(-0.1, 1.05, '(a)', transform=ax[0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')# 标注(a) (b)

# 图b绘制
ax[1].plot(x1, y1, linestyle='-', lw=lw, color="#2E5FFF", label=r'$\Gamma = 0.0005t_1$', zorder=4)
# ax[1].fill_between(x1, y1, color="#7A99FF", alpha=alpha)
ax[1].plot(x2, y2, linestyle='-', lw=lw, color="#C572FD", label=r'$\Gamma = 0.001t_1$', zorder=4)
# ax[1].fill_between(x2, y2, color="#D9A0FF", alpha=alpha)
ax[1].plot(x3, y3, linestyle='-', lw=lw, color="#FFB108", label=r'$\Gamma = 0.005t_1$', zorder=4)
# ax[1].fill_between(x3, y3, color="#FFD67A", alpha=alpha)
# ax[1].plot(x4, y4, linestyle='--', lw=lw, color="#5132FF", label=r'$\Gamma = 0.01$', zorder=4)
# ax[1].plot(x5, y5, linestyle='-.', lw=lw, color="#30ABE4", label=r'$\Gamma = 0.02$', zorder=4)
ax[1].plot(x6, y6, linestyle='-', lw=lw, color="#FF6565", label=r'$\Gamma = 0.03t_1$', zorder=4)
# ax[1].fill_between(x6, y6, color="#FF9A9A", alpha=alpha)

# 计算图b实际绘制数据的y轴最大值（在x范围[0,2]内）
x_mask_b1 = (x1 >= 0) & (x1 <= 2)
x_mask_b2 = (x2 >= 0) & (x2 <= 2)
x_mask_b3 = (x3 >= 0) & (x3 <= 2)
x_mask_b6 = (x6 >= 0) & (x6 <= 2)
y_plotted_b = np.concatenate([y1[x_mask_b1], y2[x_mask_b2], y3[x_mask_b3], y6[x_mask_b6]])
y_min_b = np.min(y_plotted_b)
y_max_b = np.max(y_plotted_b)

ax[1].set_xlabel(r'$\hbar\omega\ [t_1]$', fontsize=16)
ax[1].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$', fontsize=16)
ax[1].tick_params(axis='both', which='both', top=True, labelbottom=True, right=True, direction='in', width=0.5)
ax[1].set_xlim(0, 2)  # 设置x轴范围
ax[1].set_ylim(-0.0005, y_max_b * 1.1)
# 设置legend到左上，加粗图例线条
leg1 = ax[1].legend(handlelength=2, frameon=False, loc='upper left', fontsize=legend_fontsize)
for legline in leg1.get_lines():
    legline.set_linewidth(legend_lw)   # 设置图例线条宽度

ax[1].minorticks_on()
ax[1].tick_params(axis='both', which='major', length=3, width=0.4, direction='in') #设置主次刻度线的长度和宽度
ax[1].tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax[1].xaxis.set_tick_params(labelsize=12)
# 设置y轴刻度字体大小
ax[0].yaxis.set_tick_params(labelsize=12)
ax[1].yaxis.set_tick_params(labelsize=12) 


ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))  # x轴主刻度之间有5个次刻度
ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))  # y轴主刻度之间有5个次刻度
ax[1].text(-0.1, 1.05, '(b)', transform=ax[1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')


plt.tight_layout()  # 使用tight_layout替代subplots_adjust
plt.subplots_adjust(hspace=0)
plt.savefig('result/picture/optical conductivity 1.pdf')
plt.clf()