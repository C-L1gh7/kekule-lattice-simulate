import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from tqdm import tqdm
from numba import njit
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

import basic_function as bf

####################################################################################

# t1:inner hopping, t2:inter hopping
t1 = 2.8
t2 = 0.5*t1
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
bulkcolor = '#21467A'
edgecolor = '#DBA972'
####################################################################################

#设置存储位置   
mkpath="result/picture"# 打印文件名
bf.mkdir(mkpath)

#####################################################################################
lattice = bf.O_keku()

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

# solver2 = pb.solver.arpack(model2, k=model2.system.num_sites-1)# with corner
# eigenvalues2 = solver2.eigenvalues  # position in [nm]
# np.savetxt('matrix_output2.txt', eigenvalues2)
eigenvalues2 = np.genfromtxt('matrix_output2.txt', dtype=None)
states_num2 = list(range(len(eigenvalues2)))

corner_state = int((len(eigenvalues2)/2) - 2)
values_corner = eigenvalues2[corner_state:corner_state+6]
CornerState_index = list(range(corner_state, corner_state+6))


EdgeState_index = np.where(np.abs(eigenvalues2) <= 1)
BulkState_index = np.where(np.abs(eigenvalues2) > 1)
EdgeState_index = EdgeState_index[0].tolist()
BulkState_index = BulkState_index[0].tolist()
EdgeState = eigenvalues2[EdgeState_index]
BulkState = eigenvalues2[BulkState_index]

fig, axs = plt.subplots(1, 2)

axs[0].scatter(states_num1, eigenvalues1, s=0.1, color=bulkcolor)
axs[0].set_xlim((model1.system.num_sites/2)-1700,(model1.system.num_sites/2)+1700)
axs[0].set_ylim(-2.5,2.5)
axs[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) #隐藏x坐标轴
axs[0].tick_params(axis='y', labelsize=12)
axs[0].set_xlabel('Energy Level', fontsize=14)
axs[0].set_ylabel('Energy(eV)', fontsize=14)
axs[0].text(-0.1, 1.05, '(a)', transform=axs[0].transAxes, fontsize=12, fontweight='bold', va='top')

axs[1].scatter(EdgeState_index, EdgeState, s=0.1, color=edgecolor, label = 'Edge')
axs[1].scatter(CornerState_index, values_corner, s=0.1, color=cornercolor, label = 'Corner')
axs[1].scatter(BulkState_index, BulkState, s=0.1, color=bulkcolor, label = 'Bulk')
axs[1].set_ylim(-2.5,2.5)
axs[1].set_xlim((model2.system.num_sites/2)-1700,(model2.system.num_sites/2)+1700)
axs[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) #隐藏x坐标轴和刻度
axs[1].set_xlabel('Energy Level', fontsize=14)
axs[1].tick_params(axis='y', which='both', labelleft=False)
axs[1].text(-0.1, 1.05, '(b)', transform=axs[1].transAxes, fontsize=12, fontweight='bold', va='top')

axins = axs[1].inset_axes((0.65, 0.2, 0.3, 0.15))
axins.scatter(EdgeState_index, EdgeState, s=0.1, color=edgecolor, label = 'Edge')
axins.scatter(CornerState_index, values_corner, s=0.1, color=cornercolor, label = 'Corner')
axins.scatter(BulkState_index, BulkState, s=0.1, color=bulkcolor, label = 'Bulk')
axins.set_xlim(corner_state-1,corner_state+6)
axins.set_ylim(-0.1,0.1)
axins.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) #隐藏x坐标轴和刻度
axins.tick_params(axis='y', which='both', labelleft=False, left=False)
x1, x2 = corner_state-100, corner_state+100
y1, y2 = -5e-2, 5e-2
rect = Rectangle((x1, y1), x2-x1, y2-y1, edgecolor='black', facecolor='none', linewidth=0.7)
axs[1].add_patch(rect)  # 添加矩形框到主图

x_inset, y_inset, width_inset, height_inset= 0.65, 0.2, 0.3, 0.15
x_inset_data1, y_inset_data1 = axs[1].transData.inverted().transform(axs[1].transAxes.transform((x_inset, y_inset)))
x_inset_data2, y_inset_data2 = axs[1].transData.inverted().transform(axs[1].transAxes.transform((x_inset+width_inset, y_inset+height_inset)))

main_corners = [
    (x1, y1),
    (x2, y2)
]
inset_corners = [
    (x_inset_data1, y_inset_data1),
    (x_inset_data2, y_inset_data2)
]


# 连接主图矩形框和 inset 角点
for (x_main, y_main), (x_inset, y_inset) in zip(main_corners, inset_corners):
    line = Line2D([x_main, x_inset], [y_main, y_inset], color='black', linewidth=0.7)
    axs[1].add_line(line)  # 直接加到主图

handles, labels = axs[1].get_legend_handles_labels()
new_order = [1, 0, 2] # 更改图例显示顺序
handles = [handles[i] for i in new_order]
labels = [labels[i] for i in new_order]
axs[1].legend(handles, labels, frameon=False, markerscale=3) #显示图例，去除边框

plt.tight_layout()
# plt.show()
plt.savefig(mkpath+'/energy_band.eps')  # 保存为PNG文件
plt.clf()
def add_eigenvector(self, other):
    # 假设两个对象都有 eigenvector 属性
    self.data = np.add(self.data, other.data)
    # 其他元素保持不变，创建一个新的类实例
    return self

fig = plt.figure()
# plot wave-function of corner state
ax1 = plt.subplot(121)
probability_map1 = solver2.calc_probability(CornerState_index[0])
for i in tqdm(CornerState_index):
    if i == CornerState_index[0]:
        continue
    probability_map2 = solver2.calc_probability([i])
    probability_map1 = add_eigenvector(probability_map1,probability_map2)

probability_map1.plot(site_radius=(0.0, 0.5), cmap=[cornercolor])
ax1.set_title('')  # 隐藏标题
ax1.axis('off')

ax1.text(0.15, 0.8, '(a)', transform=ax1.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Corner', markerfacecolor=cornercolor, markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Edge', markerfacecolor=edgecolor, markersize=10)
]


ax2 = plt.subplot(122)
probability_map1 = solver2.calc_probability(EdgeState_index[0])
for i in tqdm(EdgeState_index):
    if i == EdgeState_index[0]:
        continue
    probability_map2 = solver2.calc_probability([i])
    probability_map1 = add_eigenvector(probability_map1,probability_map2)

probability_map1.plot(site_radius=(0.0, 0.5), cmap=[edgecolor])
ax2.set_title('')  # 隐藏标题

ax2.axis('off')
ax2.text(0.15, 0.8, '(b)', transform=ax2.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
plt.tight_layout()
fig.legend(handles=legend_elements, ncol=1, frameon = False, loc=(0.43, 0.25))

plt.subplots_adjust(wspace=0)
plt.savefig(mkpath+'/wavefunction.eps')
plt.clf()
  