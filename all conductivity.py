import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import AutoMinorLocator

# 物理参数设置
t1 = 2.8
t2 = 0.5 * t1
POE = 0.0
NOE = 0.0

# 晶格参数
main_L = 1.0        # 主晶格长度
bond_L = 1.0        # 键长
modulus = 2 * main_L + bond_L

# 边长参数
edgelength = 49.1   # 边长 = 4.1 + 3n (n=0,1,2,...)
edgelength1 = 48.1  # 无角情况
edgelength2 = 49.1  # 有角情况

# 物理常数
kBT = 0.01
Mu = 0
Gamma = 0.001
h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength ** 2
hS = -h_bar / S

# 文件路径设置
base_params = f"mu={Mu}kBT={kBT}gamma={Gamma}"
mkpath = f"result/anti-corner{edgelength}/{base_params}"
mkpath1 = f"result/anti-corner{edgelength1}/{base_params}"
mkpath2 = f"result/anti-corner{edgelength2}/{base_params}"

# 颜色设置
# 有角/无角颜色
with_corner_color_Between = "#FFB7F9"
without_corner_color_Between = "#FFA396"
with_corner_color = "#FF6CF3"
without_corner_color = "#FF6652"

# 不同组分颜色
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

# 绘图样式参数
alpha = 0.3
lw = 0.3
legend_fontsize = 8

# 数据加载
x = np.loadtxt(f'{mkpath}/h_omega.txt')

# 总导率数据
y_total = np.genfromtxt(f'{mkpath}/total_sigma_xx.txt', dtype=None)
y_nc = np.genfromtxt(f'{mkpath1}/total_sigma_xx.txt', dtype=None)  # 无角
y_c = np.genfromtxt(f'{mkpath2}/total_sigma_xx.txt', dtype=None)   # 有角

# 各组分导率数据
y_B_B = np.genfromtxt(f'{mkpath}/B-B_sigma_xx.txt', dtype=None)
y_B_C = np.genfromtxt(f'{mkpath}/B-C_sigma_xx.txt', dtype=None)
y_B_E = np.genfromtxt(f'{mkpath}/B-E_sigma_xx.txt', dtype=None)
y_E_E = np.genfromtxt(f'{mkpath}/E-E_sigma_xx.txt', dtype=None)
y_E_C = np.genfromtxt(f'{mkpath}/E-C_sigma_xx.txt', dtype=None)
y_C_C = np.genfromtxt(f'{mkpath}/C-C_sigma_xx.txt', dtype=None)

# 创建子图
fig, ax = plt.subplots(3, 1, sharex=True)

# 子图(a): 总导率比较
ax[0].plot(x, y_c, linestyle='--', lw=lw, color=with_corner_color, label='L=49')
ax[0].fill_between(x, y_c, color=with_corner_color_Between, alpha=alpha)
ax[0].plot(x, y_nc, linestyle=':', lw=lw, color=without_corner_color, label='L=48')
ax[0].fill_between(x, y_nc, color=without_corner_color_Between, alpha=alpha)

ax[0].set_xlim(0, 2.2)
ax[0].legend(fontsize=legend_fontsize, markerscale=1, labelspacing=0.4, 
            handlelength=2, frameon=False, loc='upper left')

# 根据xlim范围设置ylim
x_mask = (x >= 0) & (x <= 2.2)
y_c_range = y_c[x_mask]
y_nc_range = y_nc[x_mask]
y_min = min(np.min(y_c_range), np.min(y_nc_range))
y_max = max(np.max(y_c_range), np.max(y_nc_range))
y_margin = (y_max - y_min) * 0.05
ax[0].set_ylim(y_min - y_margin, y_max + y_margin)

ax[0].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')

# 子图(b): B-B和E-B组分
ax[1].plot(x, y_B_B, linestyle='--', lw=lw, color=bb, label='B-B')
ax[1].fill_between(x, y_B_B, color=bb_Between, alpha=alpha)
ax[1].plot(x, y_B_E, linestyle=':', lw=lw, color=eb, label='E-B')
ax[1].fill_between(x, y_B_E, color=eb_Between, alpha=alpha)

ax[1].set_xlim(0, 2.2)
ax[1].legend(fontsize=legend_fontsize, markerscale=1, labelspacing=0.4, 
            handlelength=2, frameon=False, loc='upper left')

# 设置ylim
y_B_B_range = y_B_B[x_mask]
y_B_E_range = y_B_E[x_mask]
y_min_1 = min(np.min(y_B_B_range), np.min(y_B_E_range))
y_max_1 = max(np.max(y_B_B_range), np.max(y_B_E_range))
y_margin_1 = (y_max_1 - y_min_1) * 0.05
ax[1].set_ylim(y_min_1 - y_margin_1, y_max_1 + y_margin_1)

ax[1].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')

# 子图(c): C-B, C-E和E-E组分
ax[2].plot(x, y_B_C, linestyle=':', lw=lw, color=cb, label='C-B')
ax[2].fill_between(x, y_B_C, color=cb_Between, alpha=alpha)
ax[2].plot(x, y_E_C, linestyle=':', lw=lw, color=ce, label='C-E')
ax[2].fill_between(x, y_E_C, color=ce_Between, alpha=alpha)
ax[2].plot(x, y_E_E, linestyle='-', lw=lw, color=ee, label='E-E')
ax[2].fill_between(x, y_E_E, color=ee_Between, alpha=alpha)

ax[2].set_xlim(0, 2.2)
ax[2].legend(fontsize=legend_fontsize, markerscale=1, labelspacing=0.4, 
            handlelength=2, frameon=False, loc='upper left')

# 设置ylim
y_B_C_range = y_B_C[x_mask]
y_E_C_range = y_E_C[x_mask]
y_E_E_range = y_E_E[x_mask]
y_min_2 = min(np.min(y_B_C_range), np.min(y_E_C_range), np.min(y_E_E_range))
y_max_2 = max(np.max(y_B_C_range), np.max(y_E_C_range), np.max(y_E_E_range))
y_margin_2 = (y_max_2 - y_min_2) * 0.05
ax[2].set_ylim(y_min_2 - y_margin_2, y_max_2 + y_margin_2)

ax[2].set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$')
ax[2].set_xlabel(r'$\hbar\omega\ [t]$')

# 设置所有子图的刻度参数
for i, a in enumerate(ax):
    # 基本刻度设置
    a.tick_params(axis='both', which='both', top=True, right=True, 
                 direction='in', width=0.5)
    a.tick_params(axis='both', which='major', length=3, width=0.4, 
                 direction='in')
    a.tick_params(axis='both', which='minor', length=2, width=0.2, 
                 direction='in')
    
    # 次刻度设置
    a.minorticks_on()
    a.xaxis.set_minor_locator(AutoMinorLocator(5))
    a.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    # 除最后一个子图外，隐藏x轴标签
    if i < len(ax) - 1:
        a.tick_params(labelbottom=False)
    else:
        # 最后一个子图显示x轴标签
        a.tick_params(labelbottom=True)
        # 设置x轴刻度，包括2.2
        a.set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.2])
        a.set_yticks([0, 0.002])

# 添加子图标签
ax[0].text(-0.05, 0.95, '(a)', transform=ax[0].transAxes, 
          fontsize=10, fontweight='bold')
ax[1].text(-0.05, 0.9, '(b)', transform=ax[1].transAxes, 
          fontsize=10, fontweight='bold')
ax[2].text(-0.05, 0.9, '(c)', transform=ax[2].transAxes, 
          fontsize=10, fontweight='bold')

# 调整布局并保存
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig('result/picture/optical conductivity 2.pdf')
# plt.show()