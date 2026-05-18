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
Gamma = 0.001
h_bar = 1.0

####################################################################################

alpha = 0.3
lw = 0.5  # 线条粗细
legend_lw = 1.2  # legend线条粗细
legend_fontsize = 14  # legend字体大小

# ==================== 图a数据路径 (L=49, edgelength=49.1) ====================
mkpath_a0   = "result/anti-corner49.1/mu=0kBT=0.01gamma=0.001"
mkpath_a02  = "result/anti-corner49.1/mu=0.2kBT=0.01gamma=0.001"
mkpath_a047 = "result/anti-corner49.1/mu=0.47kBT=0.01gamma=0.001"
mkpath_a108 = "result/anti-corner49.1/mu=1.08kBT=0.01gamma=0.001"

x_a0   = np.loadtxt(mkpath_a0   + '/h_omega.txt')
y_a0   = np.genfromtxt(mkpath_a0   + '/total_sigma_xx.txt', dtype=None)
x_a02  = np.loadtxt(mkpath_a02  + '/h_omega.txt')
y_a02  = np.genfromtxt(mkpath_a02  + '/total_sigma_xx.txt', dtype=None)
x_a047 = np.loadtxt(mkpath_a047 + '/h_omega.txt')
y_a047 = np.genfromtxt(mkpath_a047 + '/total_sigma_xx.txt', dtype=None)
x_a108 = np.loadtxt(mkpath_a108 + '/h_omega.txt')
y_a108 = np.genfromtxt(mkpath_a108 + '/total_sigma_xx.txt', dtype=None)

# ==================== 图b数据路径 (L=48, edgelength=48.1) ====================
mkpath_b0   = "result/edge48.1/t1=1_t2=2/mu=0.0_kBT=0.01_gamma=0.001"
mkpath_b02  = "result/edge48.1/t1=1_t2=2/mu=0.2_kBT=0.01_gamma=0.001"
mkpath_b047 = "result/edge48.1/t1=1_t2=2/mu=0.47_kBT=0.01_gamma=0.001"
mkpath_b108 = "result/edge48.1/t1=1_t2=2/mu=1.08_kBT=0.01_gamma=0.001"

x_b0   = np.loadtxt(mkpath_b0   + '/h_omega.txt')
y_b0   = np.genfromtxt(mkpath_b0   + '/total_sigma_xx.txt', dtype=None)
x_b02  = np.loadtxt(mkpath_b02  + '/h_omega.txt')
y_b02  = np.genfromtxt(mkpath_b02  + '/total_sigma_xx.txt', dtype=None)
x_b047 = np.loadtxt(mkpath_b047 + '/h_omega.txt')
y_b047 = np.genfromtxt(mkpath_b047 + '/total_sigma_xx.txt', dtype=None)
x_b108 = np.loadtxt(mkpath_b108 + '/h_omega.txt')
y_b108 = np.genfromtxt(mkpath_b108 + '/total_sigma_xx.txt', dtype=None)

# ==================== 绑图 ====================
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 5.6))

# 颜色与线型 (深色线条 / 浅色填充)
color_0       = "#FD3737"
color_0_fill  = "#FF6E6E"
color_02      = "#F9A710"
color_02_fill = "#FFCA25"
color_047     = "#42B47B"
color_047_fill= "#67DB9B"
color_108     = "#4A98ED"
color_108_fill= "#46B6FF"

# 图a绘制 (L=49)
ax[0].plot(x_a108, y_a108, linestyle='-', lw=lw, color=color_108, label=r'$\mu = 1.08t_1$', zorder=4)
ax[0].fill_between(x_a108, y_a108, color=color_108_fill, alpha=alpha)
ax[0].plot(x_a047, y_a047, linestyle='-', lw=lw, color=color_047, label=r'$\mu = 0.47t_1$', zorder=4)
ax[0].fill_between(x_a047, y_a047, color=color_047_fill, alpha=alpha)
ax[0].plot(x_a02,  y_a02,  linestyle='--', lw=lw, color=color_02,  label=r'$\mu = 0.2t_1$', zorder=4)
ax[0].fill_between(x_a02, y_a02, color=color_02_fill, alpha=alpha)
ax[0].plot(x_a0,   y_a0,   linestyle=':', lw=lw, color=color_0,  label=r'$\mu = 0$', zorder=4)
ax[0].fill_between(x_a0, y_a0, color=color_0_fill, alpha=alpha)

ax[0].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})/({\tilde t}^2 e^2/\hbar)$', fontsize=18)
ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False, right=True, direction='in', width=0.5)

# 设置legend到左上，加粗图例线条
handles, labels = ax[0].get_legend_handles_labels()
handles = handles[::-1]
labels = labels[::-1]
leg0 = ax[0].legend(handles, labels, handlelength=2, frameon=False,
                    loc='upper right', fontsize=legend_fontsize)
for legline in leg0.get_lines():
    legline.set_linewidth(legend_lw)

ax[0].set_xlim(0, 2.2)
ax[0].set_ylim(bottom=-0.001)
ax[0].minorticks_on()
ax[0].tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
ax[0].tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

# 图b绘制 (L=48)
ax[1].plot(x_b108, y_b108, linestyle='-', lw=lw, color=color_108, label=r'$\mu = 1.08t_1$', zorder=4)
ax[1].fill_between(x_b108, y_b108, color=color_108_fill, alpha=alpha)
ax[1].plot(x_b047, y_b047, linestyle='-', lw=lw, color=color_047, label=r'$\mu = 0.47t_1$', zorder=4)
ax[1].fill_between(x_b047, y_b047, color=color_047_fill, alpha=alpha)
ax[1].plot(x_b02,  y_b02,  linestyle='--', lw=lw, color=color_02,  label=r'$\mu = 0.2t_1$', zorder=4)
ax[1].fill_between(x_b02, y_b02, color=color_02_fill, alpha=alpha)
ax[1].plot(x_b0,   y_b0,   linestyle=':', lw=lw, color=color_0,  label=r'$\mu = 0$', zorder=4)
ax[1].fill_between(x_b0, y_b0, color=color_0_fill, alpha=alpha)

ax[1].set_xlabel(r'$\hbar\omega/t_1$', fontsize=18)
ax[1].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})/({\tilde t}^2 e^2/\hbar)$', fontsize=18)
ax[1].tick_params(axis='both', which='both', top=True, labelbottom=True, right=True, direction='in', width=0.5)
ax[1].set_xlim(0, 2.2)
ax[1].set_ylim(bottom=-0.001)

# 设置legend到左上，加粗图例线条
handles_b, labels_b = ax[1].get_legend_handles_labels()
handles_b = handles_b[::-1]
labels_b = labels_b[::-1]
leg1 = ax[1].legend(handles_b, labels_b, handlelength=2, frameon=False,
                    loc='upper right', fontsize=legend_fontsize)
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

# 隐藏图b y轴0.10的刻度标签（保留刻度线），避免遮挡(b)标注
from matplotlib.ticker import FuncFormatter
def ytick_fmt_b(val, _pos):
    if np.isclose(val, 0.10):
        return ''
    return f'{val:g}'
ax[1].yaxis.set_major_formatter(FuncFormatter(ytick_fmt_b))

# 标注(a) (b)
ax[0].text(-0.045, 0.95, '(a)', transform=ax[0].transAxes, fontsize=14, fontweight='bold')
ax[1].text(-0.045, 0.88, '(b)', transform=ax[1].transAxes, fontsize=14, fontweight='bold')

# 标注晶格尺寸
ax[0].text(0.65, 0.865, r'$L = 49\,a_0$', transform=ax[0].transAxes, fontsize=14)
ax[1].text(0.65, 0.865, r'$L = 48\,a_0$', transform=ax[1].transAxes, fontsize=14)

plt.tight_layout()
plt.subplots_adjust(hspace=0)
output_pdf = os.path.abspath('result/picture/optical_conductivity_MuCompare.pdf')
plt.savefig(output_pdf)
try:
    os.startfile(output_pdf)
except OSError as e:
    print(f"无法自动打开文件: {e}")
plt.clf()
