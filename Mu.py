import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

# ================== 参数设置 ==================
edgelength = 49.1
kBT = 0.01
Gamma = 0.001
h_bar = 1.0
alpha = 0.3
lw = 0.8
legend_lw = 1.2
legend_fontsize = 14

# ================== 数据路径 ==================
mkpath1 = f"result/anti-corner{edgelength}/mu=0.1kBT={kBT}gamma={Gamma}"
mkpath2 = f"result/anti-corner{edgelength}/mu=0.3kBT={kBT}gamma={Gamma}"
mkpath3 = f"result/anti-corner{edgelength}/mu=0.0175kBT={kBT}gamma={Gamma}"
mkpath4 = f"result/anti-corner{edgelength}/mu=0.2kBT={kBT}gamma={Gamma}"
mkpath5 = f"result/anti-corner{edgelength}/mu=1e-05kBT={kBT}gamma={Gamma}"

# ================== 读取数据 ==================
x1 = np.loadtxt(mkpath1 + '/h_omega.txt')
x2 = np.loadtxt(mkpath2 + '/h_omega.txt')
x3 = np.loadtxt(mkpath3 + '/h_omega.txt')
x4 = np.loadtxt(mkpath4 + '/h_omega.txt')   
x5 = np.loadtxt(mkpath5 + '/h_omega.txt')

y1 = np.genfromtxt(mkpath1 + '/total_sigma_xx.txt', dtype=None)
y2 = np.genfromtxt(mkpath2 + '/total_sigma_xx.txt', dtype=None)
y3 = np.genfromtxt(mkpath3 + '/total_sigma_xx.txt', dtype=None)
y4 = np.genfromtxt(mkpath4 + '/total_sigma_xx.txt', dtype=None)
y5 = np.genfromtxt(mkpath5 + '/total_sigma_xx.txt', dtype=None)

# ================== 绘制图像 ==================
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(x2, y2, '-', lw=lw, color="#FF64F7", label=r'$\mu = 0.3t_1$', zorder=4)
ax.plot(x4, y4, '--', lw=lw, color="#FF8C00", label=r'$\mu = 0.2t_1$', zorder=4)

ax.plot(x5, y5, ':', lw=lw, color="#000000", label=r'$\mu = 1e-05t_1$', zorder=4)
ax.plot(x1, y1, '-', lw=lw, color="#1DCDD3", label=r'$\mu = 0.1t_1$', zorder=4)
ax.plot(x3, y3, '-', lw=lw, color="#24863D", label=r'$\mu = 0.0175t_1$', zorder=4)



# ================== 自动计算y轴范围 ==================
# x_mask = (x1 >= 0) & (x1 <= 2.2)
# y_plotted = np.concatenate([y1[x_mask], y2[x_mask], y3[x_mask], y4[x_mask], y5[x_mask]])
# y_max = np.max(y_plotted)

# # ax.set_xlim(0, 2.2)
# ax.set_ylim(-0.0005, y_max * 1.1)
ax.set_ylabel(r'$\Re(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$', fontsize=16)

# ================== 坐标轴与刻度 ==================
ax.tick_params(axis='both', which='both', top=True, right=True,
               direction='in', width=0.5, labelsize=12)
ax.minorticks_on()
ax.tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
ax.tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

# ================== 图例与标注 ==================
handles, labels = ax.get_legend_handles_labels()
handles, labels = handles[::-1], labels[::-1]
leg = ax.legend(handles, labels, handlelength=2, frameon=False,
                loc='upper left', fontsize=legend_fontsize)
for line in leg.get_lines():
    line.set_linewidth(legend_lw)

ax.text(-0.1, 1.05, '(a)', transform=ax.transAxes, fontsize=12,
        fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.savefig('result/picture/optical_conductivity_a.pdf')
plt.close()
