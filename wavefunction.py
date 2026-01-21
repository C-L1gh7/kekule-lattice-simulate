import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from matplotlib.lines import Line2D

import basic_function as bf

# 设置全局字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

####################################################################################
# 参数设置
t1 = 1
t2 = 2*t1

main_L = 1.0
bond_L = 1.0
edgelength2 = 49.1  # edgelength = 4.1+3n(n=0,1,2,...)

####################################################################################
# 颜色
cornercolor = '#C71F2D'
edgecolor = '#DBA972'

####################################################################################
# 设置存储位置
mkpath = "result/picture"
bf.mkdir(mkpath)

####################################################################################
# 建模
lattice = bf.O_keku(t1=t1, t2=t2)

model2 = pb.Model(
    lattice,
    bf.trapezoid(r=edgelength2)
)

# 求解
solver2 = pb.solver.arpack(model2, k=model2.system.num_sites-1)
eigenvalues2 = solver2.eigenvalues
print("solver2计算完成")

# 态分类
corner_state = int((len(eigenvalues2)/2) - 2)
CornerState_index = list(range(corner_state, corner_state+6))

EdgeState_index = np.where(np.abs(eigenvalues2) <= 1)[0].tolist()

print("建模完成")

####################################################################################
# 波函数绘图

def add_eigenvector(self, other):
    self.data = np.add(self.data, other.data)
    return self

fig = plt.figure()

# plot wave-function of corner state
ax1 = plt.subplot(121)
probability_map1 = solver2.calc_probability(CornerState_index[0])
for i in tqdm(CornerState_index):
    if i == CornerState_index[0]:
        continue
    probability_map2 = solver2.calc_probability([i])
    probability_map1 = add_eigenvector(probability_map1, probability_map2)

probability_map1.plot(site_radius=(0.0, 0.7), cmap=[cornercolor])
ax1.set_title('')
ax1.axis('off')
ax1.text(0.15, 0.9, '(a)', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Corner', markerfacecolor=cornercolor, markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Edge', markerfacecolor=edgecolor, markersize=10)
]

# plot wave-function of edge state
ax2 = plt.subplot(122)
probability_map1 = solver2.calc_probability(EdgeState_index[0])
for i in tqdm(EdgeState_index):
    if i == EdgeState_index[0]:
        continue
    probability_map2 = solver2.calc_probability([i])
    probability_map1 = add_eigenvector(probability_map1, probability_map2)

probability_map1.plot(site_radius=(0.0, 0.7), cmap=[edgecolor])
ax2.set_title('')
ax2.axis('off')
ax2.text(0.15, 0.9, '(b)', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

plt.tight_layout()
fig.legend(handles=legend_elements, ncol=1, frameon=False, loc=(0.43, 0.25))
plt.subplots_adjust(wspace=0)
plt.savefig(mkpath+'/wavefunction.pdf')
plt.clf()
