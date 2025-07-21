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
mkpath="result/picture/wavefunction_test"# 打印文件名
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
solver1 = pb.solver.arpack(model1, k=model1.system.num_sites-1)# without corner
eigenvalues1 = solver1.eigenvalues  # position in [nm]
np.savetxt('matrix_output1.txt', eigenvalues1)
# eigenvalues1 = np.genfromtxt('matrix_output1.txt', dtype=None)
states_num1 = list(range(len(eigenvalues1)))

solver2 = pb.solver.arpack(model2, k=model2.system.num_sites-1)# with corner
eigenvalues2 = solver2.eigenvalues  # position in [nm]
np.savetxt('matrix_output2.txt', eigenvalues2)
# eigenvalues2 = np.genfromtxt('matrix_output2.txt', dtype=None)
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


total_samples = len(BulkState_index)
target_samples = 100
# 生成等距索引位置（在列表长度范围内等距分布的100个点）
indices = np.linspace(0, total_samples - 1, target_samples, dtype=int)
# 采样对应的值
sampled_indices = [BulkState_index[i] for i in indices]
for i in tqdm(sampled_indices):
    probability_map = solver2.calc_probability(i)
    probability_map.plot(site_radius=(0.05, 0.3))
    plt.savefig(f'{mkpath}/BulkState_{i:04d}.pdf')