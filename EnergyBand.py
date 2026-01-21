import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from numba import njit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D

import basic_function as bf


# 设置全局字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'  # 数学公式使用 STIX 字体(与 Times New Roman 风格一致)

####################################################################################

# t1:inner hopping, t2:inter hopping
t1 = 1
t2 = 2*t1
POE = 0.0
NOE = 0.0

main_L = 1.0  # length of main lattice
bond_L = 1.0  # bond length between lattices.
modulus = 2 * main_L + bond_L
edgelength1 = 48.1 # edgelength = 4.1+3n(n=0,1,2,...)
edgelength2 = 49.1 # edgelength = 4.1+3n(n=0,1,2,...)
####################################################################################

kBT = 0.01
Mu = 0
Gamma = 0.001
h_bar = 1.0

####################################################################################
# color
cornercolor = '#C71F2D'
bulkcolor = "#38557E"
edgecolor = '#DBA972'

s = 0.7  # marker size for scatter plot
####################################################################################

#设置存储位置   
mkpath="result/picture"# 打印文件名
bf.mkdir(mkpath)

#####################################################################################
lattice = bf.O_keku(t1=t1, t2=t2)

model1 = pb.Model(
    lattice,
    bf.trapezoid(r=edgelength1)
)

model2 = pb.Model(
    lattice,
    bf.trapezoid(r=edgelength2)
)

# 绘制energyband
# solver1 = pb.solver.arpack(model1, k=model1.system.num_sites-1)# without corner
# eigenvalues1 = solver1.eigenvalues  # position in [nm]
# np.savetxt('matrix_output1.txt', eigenvalues1)
eigenvalues1 = np.genfromtxt('matrix_output1.txt', dtype=None)
states_num1 = list(range(len(eigenvalues1)))
print("solver1计算完成")

# solver2 = pb.solver.arpack(model2, k=model2.system.num_sites-1)# with corner
# eigenvalues2 = solver2.eigenvalues  # position in [nm]
# np.savetxt('matrix_output2.txt', eigenvalues2)
eigenvalues2 = np.genfromtxt('matrix_output2.txt', dtype=None)
states_num2 = list(range(len(eigenvalues2)))
print("solver2计算完成")
corner_state = int((len(eigenvalues2)/2) - 2)
values_corner = eigenvalues2[corner_state:corner_state+6]
CornerState_index = list(range(corner_state, corner_state+6))


EdgeState_index = np.where(np.abs(eigenvalues2) <= 1)
BulkState_index = np.where(np.abs(eigenvalues2) > 1)
EdgeState_index = EdgeState_index[0].tolist()
BulkState_index = BulkState_index[0].tolist()
EdgeState = eigenvalues2[EdgeState_index]
BulkState = eigenvalues2[BulkState_index]

print("建模完成")

fig, axs = plt.subplots(1, 2)

# 原来的图(b)现在变成图(a) - 左边的图
axs[0].scatter(EdgeState_index, EdgeState, s=s, color=edgecolor, label = 'Edge')
axs[0].scatter(CornerState_index, values_corner, s=s, color=cornercolor, label = 'Corner')
axs[0].scatter(BulkState_index, BulkState, s=s, color=bulkcolor, label = 'Bulk')
axs[0].set_ylim(-2.5,2.5)
axs[0].set_xlim((model2.system.num_sites/2)-1700,(model2.system.num_sites/2)+1700)
axs[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) #隐藏x坐标轴和刻度
axs[0].tick_params(axis='y', labelsize=12, direction='in')  # 刻度向内
axs[0].set_xlabel('Energy Level', fontsize=14)
axs[0].set_ylabel(r'$\mathrm{Energy}[t_1]$', fontsize=14)
axs[0].text(-0.1, 1.05, '(a)', transform=axs[0].transAxes, fontsize=12, fontweight='bold', va='top')

# 为图(a)添加inset图
axins = axs[0].inset_axes((0.65, 0.2, 0.3, 0.15))
axins.scatter(EdgeState_index, EdgeState, s=s, color=edgecolor, label = 'Edge')
axins.scatter(CornerState_index, values_corner, s=s, color=cornercolor, label = 'Corner')
axins.scatter(BulkState_index, BulkState, s=s, color=bulkcolor, label = 'Bulk')
axins.set_xlim(corner_state-1,corner_state+6)
axins.set_ylim(-0.0015,0.0015)
axins.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))#使用科学计数法显示y轴刻度
axins.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) #隐藏x坐标轴和刻度
axins.tick_params(axis='y', labelsize=8, direction='in')  # 添加纵坐标刻度，刻度向内
axins.yaxis.set_major_locator(plt.MaxNLocator(3))  # 设置纵坐标刻度数量
x1, x2 = corner_state-100, corner_state+100
y1, y2 = -5e-2, 5e-2
# 主图中的虚线框
rect = Rectangle((x1, y1), x2-x1, y2-y1, edgecolor='red', facecolor='none', 
                 linewidth=0.7, linestyle='--')
axs[0].add_patch(rect)

# ==================== 添加带间隙的箭头 ====================
# 虚线框右下角（数据坐标）
p1_data = (x2, y1)

# 需要先绘制图形以获取正确的坐标转换
fig.canvas.draw()

# 将虚线框右下角从数据坐标转换为 axes 坐标
p1_display = axs[0].transData.transform(p1_data)
p1_axes = axs[0].transAxes.inverted().transform(p1_display)

# inset 左上角在 axes 坐标系中的位置
# inset_axes 的位置是 (0.65, 0.2, 0.3, 0.15)，左上角 y = 0.2 + 0.15
p2_axes = (0.65, 0.35)

# 计算方向向量
dx = p2_axes[0] - p1_axes[0]
dy = p2_axes[1] - p1_axes[1]
length = np.sqrt(dx**2 + dy**2)
dx_norm = dx / length
dy_norm = dy / length

# 设置间隙（绝对距离，在 axes 坐标系中）
gap_start = 0.03  # 起点处的间隙
gap_end = 0.02    # 终点处的间隙

# 计算缩短后的起点和终点（在 axes 坐标系中）
arrow_start = (p1_axes[0] + gap_start * dx_norm, p1_axes[1] + gap_start * dy_norm)
arrow_end = (p2_axes[0] - gap_end * dx_norm, p2_axes[1] - gap_end * dy_norm)

# 绘制箭头
arrow = FancyArrowPatch(
    arrow_start, arrow_end,
    transform=axs[0].transAxes,
    arrowstyle='->,head_length=3,head_width=2',
    color='red',
    linewidth=0.9,
    linestyle='-'
)
axs[0].add_patch(arrow)
# ==================== 箭头添加完成 ====================

handles, labels = axs[0].get_legend_handles_labels()
new_order = [1, 0, 2] # 更改图例显示顺序
handles = [handles[i] for i in new_order]
labels = [labels[i] for i in new_order]
axs[0].legend(handles, labels, frameon=False, markerscale=3) #显示图例，去除边框

# 原来的图(a)现在变成图(b) - 右边的图
axs[1].scatter(states_num1, eigenvalues1, s=0.1, color=bulkcolor)
axs[1].set_xlim((model1.system.num_sites/2)-1700,(model1.system.num_sites/2)+1700)
axs[1].set_ylim(-2.5,2.5)
axs[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) #隐藏x坐标轴
axs[1].tick_params(axis='y', which='both', labelleft=False, direction='in')  # 刻度向内
axs[1].set_xlabel('Energy Level', fontsize=14)
axs[1].text(-0.1, 1.05, '(b)', transform=axs[1].transAxes, fontsize=12, fontweight='bold', va='top')

plt.tight_layout()
# plt.show()
plt.savefig(mkpath+'/energy_band.pdf')
  