import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from tqdm import tqdm
from matplotlib.ticker import AutoMinorLocator

import basic_function as bf

# 设置全局字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'  # 数学公式使用 STIX 字体(与 Times New Roman 风格一致)

####################################################################################
kBT = 0.01
Mu = 0
h_bar = 1.0

####################################################################################

lw = 0.5  # 线条粗细
legend_lw = 1.2  # legend线条粗细
legend_fontsize = 14  # legend字体大小

# 颜色定义
color_0005 = "#2E5FFF"
color_001  = "#C572FD"
color_005  = "#FFB108"
color_03   = "#FF6565"

# ==================== 图a数据路径 (L=49, edgelength=49.1) ====================
mkpath_a1 = "result/anti-corner49.1/mu=0kBT=0.01gamma=0.0005"
mkpath_a2 = "result/anti-corner49.1/mu=0kBT=0.01gamma=0.001"
mkpath_a3 = "result/anti-corner49.1/mu=0kBT=0.01gamma=0.005"
mkpath_a6 = "result/anti-corner49.1/mu=0kBT=0.01gamma=0.03"

x_a1 = np.loadtxt(mkpath_a1 + '/h_omega.txt')
y_a1 = np.genfromtxt(mkpath_a1 + '/total_sigma_xx.txt', dtype=None)
x_a2 = np.loadtxt(mkpath_a2 + '/h_omega.txt')
y_a2 = np.genfromtxt(mkpath_a2 + '/total_sigma_xx.txt', dtype=None)
x_a3 = np.loadtxt(mkpath_a3 + '/h_omega.txt')
y_a3 = np.genfromtxt(mkpath_a3 + '/total_sigma_xx.txt', dtype=None)
x_a6 = np.loadtxt(mkpath_a6 + '/h_omega.txt')
y_a6 = np.genfromtxt(mkpath_a6 + '/total_sigma_xx.txt', dtype=None)

# ==================== 图b数据路径 (L=48, edgelength=48.1) ====================
mkpath_b1 = "result/edge48.1/t1=1_t2=2/mu=0_kBT=0.01_gamma=0.0005"
mkpath_b2 = "result/edge48.1/t1=1_t2=2/mu=0_kBT=0.01_gamma=0.001"
mkpath_b3 = "result/edge48.1/t1=1_t2=2/mu=0_kBT=0.01_gamma=0.005"
mkpath_b6 = "result/edge48.1/t1=1_t2=2/mu=0_kBT=0.01_gamma=0.03"

x_b1 = np.loadtxt(mkpath_b1 + '/h_omega.txt')
y_b1 = np.genfromtxt(mkpath_b1 + '/total_sigma_xx.txt', dtype=None)
x_b2 = np.loadtxt(mkpath_b2 + '/h_omega.txt')
y_b2 = np.genfromtxt(mkpath_b2 + '/total_sigma_xx.txt', dtype=None)
x_b3 = np.loadtxt(mkpath_b3 + '/h_omega.txt')
y_b3 = np.genfromtxt(mkpath_b3 + '/total_sigma_xx.txt', dtype=None)
x_b6 = np.loadtxt(mkpath_b6 + '/h_omega.txt')
y_b6 = np.genfromtxt(mkpath_b6 + '/total_sigma_xx.txt', dtype=None)

# ==================== 绑图 ====================
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 5.12))

# 图a绘制 (L=49)
ax[0].plot(x_a1, y_a1, linestyle='-', lw=lw, color=color_0005, label=r'$\Gamma = 0.0005t_1$', zorder=4)
ax[0].plot(x_a2, y_a2, linestyle='-', lw=lw, color=color_001, label=r'$\Gamma = 0.001t_1$', zorder=4)
ax[0].plot(x_a3, y_a3, linestyle='-', lw=lw, color=color_005, label=r'$\Gamma = 0.005t_1$', zorder=4)
ax[0].plot(x_a6, y_a6, linestyle='-', lw=lw, color=color_03, label=r'$\Gamma = 0.03t_1$', zorder=4)

# 计算x∈[0,2.2]范围内y轴最大值
y_max_a = max(
    np.max(y_a1[(x_a1 >= 0) & (x_a1 <= 2.2)]),
    np.max(y_a2[(x_a2 >= 0) & (x_a2 <= 2.2)]),
    np.max(y_a3[(x_a3 >= 0) & (x_a3 <= 2.2)]),
    np.max(y_a6[(x_a6 >= 0) & (x_a6 <= 2.2)]),
)

ax[0].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})/({\tilde t}^2 e^2/\hbar)$', fontsize=18)
ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False, right=True, direction='in', width=0.5)

# 设置legend，加粗图例线条
leg0 = ax[0].legend(handlelength=2, frameon=False, loc='upper left', fontsize=legend_fontsize)
for legline in leg0.get_lines():
    legline.set_linewidth(legend_lw)

ax[0].set_xlim(0, 2.2)
ax[0].set_ylim(-0.001, y_max_a * 1.1)
ax[0].minorticks_on()
ax[0].tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
ax[0].tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

# 图b绘制 (L=48)
ax[1].plot(x_b1, y_b1, linestyle='-', lw=lw, color=color_0005, label=r'$\Gamma = 0.0005t_1$', zorder=4)
ax[1].plot(x_b2, y_b2, linestyle='-', lw=lw, color=color_001, label=r'$\Gamma = 0.001t_1$', zorder=4)
ax[1].plot(x_b3, y_b3, linestyle='-', lw=lw, color=color_005, label=r'$\Gamma = 0.005t_1$', zorder=4)
ax[1].plot(x_b6, y_b6, linestyle='-', lw=lw, color=color_03, label=r'$\Gamma = 0.03t_1$', zorder=4)

ax[1].set_xlabel(r'$\hbar\omega/t_1$', fontsize=18)
ax[1].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})/({\tilde t}^2 e^2/\hbar)$', fontsize=18)
ax[1].tick_params(axis='both', which='both', top=True, labelbottom=True, right=True, direction='in', width=0.5)
ax[1].set_xlim(0, 2.2)
ax[1].set_ylim(-0.001, y_max_a * 1.1)

# 设置legend，加粗图例线条
leg1 = ax[1].legend(handlelength=2, frameon=False, loc='upper left', fontsize=legend_fontsize)
for legline in leg1.get_lines():
    legline.set_linewidth(legend_lw)

ax[1].minorticks_on()
ax[1].tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
ax[1].tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax[1].xaxis.set_tick_params(labelsize=14)
ax[0].yaxis.set_tick_params(labelsize=14)
ax[1].yaxis.set_tick_params(labelsize=14)

ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))

# 在x轴添加2.2的tick
xticks = [t for t in ax[1].get_xticks() if t <= 2.2]
if 2.2 not in xticks:
    xticks.append(2.2)
ax[1].set_xticks(xticks)
ax[1].set_xlim(0, 2.2)

# 标注(a) (b)
ax[0].text(-0.045, 0.95, '(a)', transform=ax[0].transAxes, fontsize=14, fontweight='bold')
ax[1].text(-0.045, 0.922, '(b)', transform=ax[1].transAxes, fontsize=14, fontweight='bold')

# 标注晶格尺寸（右上角）
ax[0].text(0.97, 0.88, r'$L = 49\,a_0$', transform=ax[0].transAxes, fontsize=14, ha='right')
ax[1].text(0.97, 0.88, r'$L = 48\,a_0$', transform=ax[1].transAxes, fontsize=14, ha='right')

plt.tight_layout()
plt.subplots_adjust(hspace=0)
output_pdf = os.path.abspath('result/picture/optical conductivity GammaCompare.pdf')
plt.savefig(output_pdf)
try:
    os.startfile(output_pdf)
except OSError as e:
    print(f"无法自动打开文件: {e}")
plt.clf()
