import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator

# 设置全局字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'  # 数学公式使用 STIX 字体(与 Times New Roman 风格一致)


def decimate_data(x, y, max_points=5000):
    """对数据进行降采样，保留曲线特征"""
    if len(x) <= max_points:
        return x, y
    step = len(x) // max_points
    if step < 1:
        step = 1
    indices = np.arange(0, len(x), step)
    return x[indices], y[indices]

# ==================== 参数设置区域 ====================
edgelength1 = 48.1    # L=48 (无角态)
edgelength2 = 49.1    # L=49 (有角态)
Mu = 1.08
kBT = 0.01
Gamma = 0.001
t1 = 1
t2 = 2
# ====================================================

# 路径
mkpath = f"result/anti-corner{edgelength2}/mu={Mu}kBT={kBT}gamma={Gamma}"
mkpath1 = f"result/edge{edgelength1}/t1={t1}_t2={t2}/mu={Mu}_kBT={kBT}_gamma={Gamma}"
mkpath2 = f"result/anti-corner{edgelength2}/mu={Mu}kBT={kBT}gamma={Gamma}"

# 颜色
with_corner_color_Between = "#FFB7F9"
without_corner_color_Between = "#FFA396"
with_corner_color = "#FF6CF3"
without_corner_color = "#FF6652"

bb_Between = "#8AC7FF"; cb_Between = "#F2DBFD"; eb_Between = "#92FFF0"
ee_Between = "#9DFFBA"; ce_Between = "#FFFD72"; cc_Between = "#FF9C9C"
bb = "#40A3FF"; cb = "#CC5FFF"; eb = "#24D8C0"
ee = "#1BA846"; ce = "#D6C800"; cc = "#EC5252"

alpha = 0.3
lw = 0.3
legend_fontsize = 10  # legend字体大小
legend_lw = 1  # legend线条粗细


# 数据
x = np.loadtxt(mkpath + '/h_omega.txt')
y_total = np.genfromtxt(mkpath + '/total_sigma_xx.txt', dtype=None)
y_B_B = np.genfromtxt(mkpath + '/B-B_sigma_xx.txt', dtype=None)
y_B_C = np.genfromtxt(mkpath + '/B-C_sigma_xx.txt', dtype=None)
y_B_E = np.genfromtxt(mkpath + '/B-E_sigma_xx.txt', dtype=None)
y_E_E = np.genfromtxt(mkpath + '/E-E_sigma_xx.txt', dtype=None)
y_E_C = np.genfromtxt(mkpath + '/E-C_sigma_xx.txt', dtype=None)
y_C_C = np.genfromtxt(mkpath + '/C-C_sigma_xx.txt', dtype=None)
x_nc = np.genfromtxt(mkpath1 + '/h_omega.txt', dtype=None)
y_nc = np.genfromtxt(mkpath1 + '/total_sigma_xx.txt', dtype=None)
x_c = np.genfromtxt(mkpath2 + '/h_omega.txt', dtype=None)
y_c = np.genfromtxt(mkpath2 + '/total_sigma_xx.txt', dtype=None)

# 数据降采样（减少数据点以提高绘图性能）
# 重要：每组数据需要独立降采样，避免x轴数据混乱
x_orig = x.copy()  # 保存原始x
x, y_total = decimate_data(x_orig, y_total)
x_bb, y_B_B = decimate_data(x_orig, y_B_B)
x_bc, y_B_C = decimate_data(x_orig, y_B_C)
x_be, y_B_E = decimate_data(x_orig, y_B_E)
x_ee, y_E_E = decimate_data(x_orig, y_E_E)
x_ec, y_E_C = decimate_data(x_orig, y_E_C)
x_cc, y_C_C = decimate_data(x_orig, y_C_C)
x_nc, y_nc = decimate_data(x_nc, y_nc)
x_c, y_c = decimate_data(x_c, y_c)

# 三个子图
fig, ax = plt.subplots(3, 1, sharex=True)

# ===== 图 (a)：L=49 vs L=48 总光电导 =====
ax[0].plot(x_c, y_c, linestyle='--', lw=lw, color=with_corner_color, label='L=49')
ax[0].fill_between(x_c, y_c, color=with_corner_color_Between, alpha=alpha)
ax[0].plot(x_nc, y_nc, linestyle=':', lw=lw, color=without_corner_color, label='L=48')
ax[0].fill_between(x_nc, y_nc, color=without_corner_color_Between, alpha=alpha)
leg0 = ax[0].legend(handlelength=2, frameon=False, loc='upper left', fontsize=legend_fontsize)
for legline in leg0.get_lines():
    legline.set_linewidth(legend_lw)
ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False, right=True, direction='in', width=0.5)
ax[0].set_xlim(0, 2.2)

# ===== 图 (b)：C-B, C-E, E-E =====
ax[1].plot(x_bc, y_B_C, linestyle=':', lw=lw, color=cb, label='C-B')
ax[1].fill_between(x_bc, y_B_C, color=cb_Between, alpha=alpha)
ax[1].plot(x_ec, y_E_C, linestyle=':', lw=lw, color=ce, label='C-E')
ax[1].fill_between(x_ec, y_E_C, color=ce_Between, alpha=alpha)
ax[1].plot(x_ee, y_E_E, linestyle='-', lw=lw, color=ee, label='E-E')
ax[1].fill_between(x_ee, y_E_E, color=ee_Between, alpha=alpha)
leg1 = ax[1].legend(fontsize=legend_fontsize, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
for legline in leg1.get_lines():
    legline.set_linewidth(legend_lw)
ax[1].minorticks_on()
ax[1].set_xlim(0, 2.2)
ax[1].tick_params(axis='both', which='both', top=True, labelbottom=False, right=True, direction='in', width=0.5)
ax[1].set_yticks([0.000, 0.002])
ax[1].set_yticklabels(['0.000', '0.002'])
ax[1].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})/({\tilde t}^2 e^2/\hbar)$', fontsize=14)


# ===== 图 (c)：B-B, E-B =====
ax[2].plot(x_bb, y_B_B, linestyle='--', lw=lw, color=bb, label='B-B')
ax[2].fill_between(x_bb, y_B_B, color=bb_Between, alpha=alpha)
ax[2].plot(x_be, y_B_E, linestyle=':', lw=lw, color=eb, label='E-B')
ax[2].fill_between(x_be, y_B_E, color=eb_Between, alpha=alpha)
leg2 = ax[2].legend(fontsize=legend_fontsize, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
for legline in leg2.get_lines():
    legline.set_linewidth(legend_lw)
ax[2].minorticks_on()
ax[2].set_xlim(0, 2.2)
ax[2].set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.2])
ax[2].tick_params(axis='both', which='both', top=True, labelbottom=True, right=True, direction='in', width=0.5)
ax[2].set_xlabel(r'$\hbar\omega/t_1$', fontsize=14)

ax[2].xaxis.set_tick_params(labelsize=12)

# 设置y轴刻度字体大小
ax[0].yaxis.set_tick_params(labelsize=12)
ax[1].yaxis.set_tick_params(labelsize=12)
ax[2].yaxis.set_tick_params(labelsize=12)

# 统一刻度样式
for a in ax:
    a.tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
    a.tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
    a.xaxis.set_minor_locator(AutoMinorLocator(5))
    a.yaxis.set_minor_locator(AutoMinorLocator(5))

# 标注 (a)(b)(c)
ax[0].text(-0.065, 0.9, '(a)', transform=ax[0].transAxes, fontsize=14, fontweight='bold')
ax[1].text(-0.065, 0.87, '(b)', transform=ax[1].transAxes, fontsize=14, fontweight='bold')
ax[2].text(-0.065, 0.87, '(c)', transform=ax[2].transAxes, fontsize=14, fontweight='bold')

# 调整排版
plt.tight_layout()
plt.subplots_adjust(hspace=0)

output_pdf = os.path.abspath('result/picture/optical_conductivity_mu.pdf')
plt.savefig(output_pdf)
try:
    os.startfile(output_pdf)
except OSError as e:
    print(f"无法自动打开文件: {e}")
