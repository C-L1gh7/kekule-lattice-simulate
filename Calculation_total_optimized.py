# Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import math
from tqdm import tqdm

import basic_function as bf

####################################################################################

# t1:inner hopping, t2:inter hopping
t1 = 1
t2 = 2*t1
POE = 0.0
NOE = 0.0

main_L = 1.0  # length of main lattice
bond_L = 1.0  # bond length between lattices
modulus = 2 * main_L + bond_L
edgelength = 48.1 # edgelength = 4.1+3n(n=0,1,2,...)
####################################################################################

kBT = 0.01
Mu = 0
Gamma = 0.001
h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength ** 2.
hS = -h_bar / S

####################################################################################

# 设置存储位置
mkpath = "result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)
# mkpath = "result/anti-corner"+str(edgelength)+"_"+str(t2/t1)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)
bf.mkdir(mkpath)

####################################################################################
lattice = bf.O_keku(t1=t1, t2=t2)
lattice.plot()
plt.savefig(mkpath+'/lattice.eps')
plt.close()

model = pb.Model(
    lattice,
    bf.trapezoid(r=edgelength)
)

model.plot()
plt.savefig(mkpath+'/model.eps')
plt.close()
bf.probability(model, mkpath)
plt.close('all')  # 关闭 probability 函数内产生的所有空白 figure

# 计算速度算符
length_x = len(model.system.x)
r_x = np.diag(model.system.x).astype(np.float64)  # 直接用 np.diag 创建对角矩阵
H = model.hamiltonian.todense()
V_x = (r_x @ H) - (H @ r_x)
np.save(mkpath+'/Vx.npy', V_x)  # 用二进制格式保存，更快

r_y = np.diag(model.system.y).astype(np.float64)
V_y = (r_y @ H) - (H @ r_y)
np.save(mkpath+'/Vy.npy', V_y)

solver = pb.solver.lapack(model)
values = solver.eigenvalues
vectors = solver.eigenvectors.astype(np.float64)
####################################################################################

indices_to_keep = np.where(np.abs(values) <= 2.5)  # 计算的能量范围
values = values[indices_to_keep]
vectors = vectors[indices_to_keep]  # 本征向量按行存储，与原版一致

####################################################################################
# 预计算常量，避免循环内重复计算

# 预计算费米分布
fF = 1.0 / (np.exp((values - Mu) / kBT) + 1.0)

# 预计算速度矩阵元 <n|V_x|k>，形状为 (num_states, num_states)
# 原版: aA = vectors[n].T @ V_x, inner_product_1 = aA @ vectors[k]
# 即 <n|V_x|k> = vectors[n] @ V_x @ vectors[k].T
V_x_matrix = vectors @ V_x @ vectors.T

# 预计算 |<n|V_x|k>|^2 = <n|V_x|k> * <k|V_x|n>
V_x_squared = np.array(V_x_matrix) * np.array(V_x_matrix.T)

# 预计算能量差矩阵 epsilon[n,k] = E_n - E_k
E_n, E_k = np.meshgrid(values, values, indexing='ij')
epsilon = E_n - E_k

# 预计算费米差矩阵 delta_fF[n,k] = fF_n - fF_k
fF_n, fF_k = np.meshgrid(fF, fF, indexing='ij')
delta_fF = fF_n - fF_k

# 创建掩码排除对角元 (n == k)
mask = ~np.eye(len(values), dtype=bool)

####################################################################################

def conductivity_vectorized(h_omega):
    """向量化计算光电导 - 比循环快 10-100 倍"""
    # 计算分母: epsilon * (h_omega + epsilon) + i*Gamma
    denominator = epsilon * (h_omega + epsilon) + 1j * Gamma

    # 计算 sigma，只对非对角元求和
    sigma_matrix = delta_fF * V_x_squared / denominator
    sigma = np.sum(sigma_matrix[mask])

    return sigma * 1j * hS


####################################################################################

# 绘图和状态分类（保持原有逻辑）
corner_state = int((len(values)/2) - 3)
values_corner = values[corner_state:corner_state+6]
CornerState_index = list(range(corner_state, corner_state+6))

EdgeState_index = np.where(np.abs(values) <= 1)[0].tolist()
BulkState_index = np.where(np.abs(values) > 1)[0].tolist()

plt.figure()
plt.scatter(CornerState_index, values[CornerState_index])
plt.scatter(EdgeState_index, values[EdgeState_index])
plt.scatter(BulkState_index, values[BulkState_index])
plt.savefig(mkpath+'/states_classification.png')
plt.close()

####################################################################################

# 计算光电导
# 基础采样：0-2.5005，步长 0.0005
omega_base = np.arange(0, 2.5005, 0.0005)

# # 峰位置配置（与 all conductivity.py 中的 peak_config 一致）
# peak_positions = [
#     0.39428 + 0.002,   # C-E
#     0.788596 + 0.002,  # E-E
#     1.00394 + 0.002,   # C-B
#     1.39823 + 0.002,   # E-B
#     2.00789 + 0.002    # B-B
# ]
#
# # 在每个峰附近进行密集采样（步长 0.0001）
# dense_range = 0.005  # 在峰位置左右各 0.005 范围内密集采样
# dense_step = 0.0001  # 密集采样步长
#
# x_dense = []
# for peak_pos in peak_positions:
#     # 在每个峰附近添加密集采样点
#     x_dense.append(np.arange(peak_pos - dense_range, peak_pos + dense_range, dense_step))
#
# # 合并基础采样和密集采样，去重并排序
# omega_range = np.unique(np.sort(np.concatenate([omega_base] + x_dense)))

omega_range = omega_base

np.save(mkpath+'/h_omega.npy', omega_range)
np.savetxt(mkpath+'/h_omega.txt', omega_range)

y_total = []
for omega in tqdm(omega_range, desc="计算光电导"):
    sigma = conductivity_vectorized(omega)
    y_total.append(sigma)

y_total = np.array(y_total)

# 储存数据
np.save(mkpath+'/total_sigma_xx.npy', y_total)

# 同时保存txt格式方便查看
np.savetxt(mkpath+'/total_sigma_xx.txt', y_total)

# 绘制光电导图
from matplotlib.ticker import AutoMinorLocator

# 设置全局字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

# 绘图参数
line_color = "#FF6CF3"
fill_color = "#FFB7F9"
alpha = 0.3
lw = 0.8
legend_fontsize = 10

fig, ax = plt.subplots(figsize=(8, 5))

# 绘制实部
y_real = np.real(y_total)
ax.plot(omega_range, y_real, linestyle='-', lw=lw, color=line_color, label=r'$\mathrm{Re}(\sigma_{xx})$')
ax.fill_between(omega_range, y_real, color=fill_color, alpha=alpha)

# 设置坐标轴
ax.set_xlim(0, omega_range[-1])
ax.set_xlabel(r'$\hbar\omega\ [t_1]$', fontsize=14)
ax.set_ylabel(r'$\mathrm{Re}(\sigma_{xx})\ [{\tilde t}^2 e^2/\hbar]$', fontsize=14)

# 设置刻度
ax.tick_params(axis='both', which='both', top=True, right=True, direction='in', width=0.5)
ax.tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
ax.tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

# 图例
leg = ax.legend(fontsize=legend_fontsize, frameon=False, loc='upper right')

plt.tight_layout()
plt.savefig(mkpath+'/sigma_xx.pdf')
plt.show()

print("计算完成！")
