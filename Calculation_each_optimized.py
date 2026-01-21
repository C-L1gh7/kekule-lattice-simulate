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
edgelength = 49.1 # edgelength = 4.1+3n(n=0,1,2,...)
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
r_x = np.diag(model.system.x).astype(np.float64)
H = model.hamiltonian.todense()
V_x = (r_x @ H) - (H @ r_x)
np.save(mkpath+'/Vx.npy', V_x)

r_y = np.diag(model.system.y).astype(np.float64)
V_y = (r_y @ H) - (H @ r_y)
np.save(mkpath+'/Vy.npy', V_y)

solver = pb.solver.lapack(model)
values = solver.eigenvalues
vectors = solver.eigenvectors.astype(np.float64)
####################################################################################

indices_to_keep = np.where(np.abs(values) <= 2.5)
values = values[indices_to_keep]
vectors = vectors[indices_to_keep]

####################################################################################
# 状态分类
corner_state = int((len(values)/2) - 3)
CornerState_index = list(range(corner_state, corner_state+6))
EdgeState_index = np.where(np.abs(values) <= 1)[0].tolist()
BulkState_index = np.where(np.abs(values) > 1)[0].tolist()

# 保存状态分类图
plt.figure()
plt.scatter(CornerState_index, values[CornerState_index])
plt.scatter(EdgeState_index, values[EdgeState_index])
plt.scatter(BulkState_index, values[BulkState_index])
plt.savefig(mkpath+'/states_classification.png')
plt.close()

####################################################################################
# 预计算常量

# 费米分布
fF = 1.0 / (np.exp((values - Mu) / kBT) + 1.0)

# 速度矩阵元 <n|V_x|k>
V_x_matrix = vectors @ V_x @ vectors.T

# |<n|V_x|k>|^2
V_x_squared = np.array(V_x_matrix) * np.array(V_x_matrix.T)

# 能量差矩阵
E_n, E_k = np.meshgrid(values, values, indexing='ij')
epsilon = E_n - E_k

# 费米差矩阵
fF_n, fF_k = np.meshgrid(fF, fF, indexing='ij')
delta_fF = fF_n - fF_k

# 排除对角元的掩码
mask_diag = ~np.eye(len(values), dtype=bool)

####################################################################################
# 创建状态分类掩码矩阵

n_states = len(values)

# 创建布尔掩码：标记每个态属于哪类
is_corner = np.zeros(n_states, dtype=bool)
is_edge = np.zeros(n_states, dtype=bool)
is_bulk = np.zeros(n_states, dtype=bool)

is_corner[CornerState_index] = True
is_edge[EdgeState_index] = True
is_bulk[BulkState_index] = True

# 注意：Corner态是Edge态的子集，需要区分纯Edge态
is_pure_edge = is_edge & ~is_corner

# 创建跃迁掩码矩阵 (n_states × n_states)
# mask_BB[n, k] = True 表示 n 和 k 都是 Bulk 态
mask_BB = np.outer(is_bulk, is_bulk)
mask_CC = np.outer(is_corner, is_corner)
mask_BC = (np.outer(is_bulk, is_corner) | np.outer(is_corner, is_bulk))
mask_BE = (np.outer(is_bulk, is_pure_edge) | np.outer(is_pure_edge, is_bulk))
mask_EE = np.outer(is_pure_edge, is_pure_edge)
mask_EC = (np.outer(is_pure_edge, is_corner) | np.outer(is_corner, is_pure_edge))

# 排除对角元
mask_BB &= mask_diag
mask_CC &= mask_diag
mask_BC &= mask_diag
mask_BE &= mask_diag
mask_EE &= mask_diag
mask_EC &= mask_diag

####################################################################################

def conductivity_vectorized(h_omega):
    """向量化计算各态贡献的光电导"""
    # 计算分母
    denominator = epsilon * (h_omega + epsilon) + 1j * Gamma

    # 计算基础 sigma 矩阵
    sigma_matrix = delta_fF * V_x_squared / denominator

    # 对不同态的跃迁分别求和
    sigma_total = np.sum(sigma_matrix[mask_diag])
    sigma_BB = np.sum(sigma_matrix[mask_BB])
    sigma_CC = np.sum(sigma_matrix[mask_CC])
    sigma_BC = np.sum(sigma_matrix[mask_BC])
    sigma_BE = np.sum(sigma_matrix[mask_BE])
    sigma_EE = np.sum(sigma_matrix[mask_EE])
    sigma_EC = np.sum(sigma_matrix[mask_EC])

    # 乘以系数
    factor = 1j * hS
    return (sigma_total * factor, sigma_BB * factor, sigma_BC * factor,
            sigma_BE * factor, sigma_EE * factor, sigma_EC * factor, sigma_CC * factor)

####################################################################################

# 基础采样：0-2.5005，步长 0.0005
omega_base = np.arange(0, 2.20005, 0.00005)

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

# 计算光电导
y_total = []
y_B_B = []
y_B_C = []
y_B_E = []
y_E_E = []
y_E_C = []
y_C_C = []

for omega in tqdm(omega_range, desc="计算各态光电导"):
    sigma_vals = conductivity_vectorized(omega)
    y_total.append(sigma_vals[0])
    y_B_B.append(sigma_vals[1])
    y_B_C.append(sigma_vals[2])
    y_B_E.append(sigma_vals[3])
    y_E_E.append(sigma_vals[4])
    y_E_C.append(sigma_vals[5])
    y_C_C.append(sigma_vals[6])

# 转为数组
y_total = np.array(y_total)
y_B_B = np.array(y_B_B)
y_B_C = np.array(y_B_C)
y_B_E = np.array(y_B_E)
y_E_E = np.array(y_E_E)
y_E_C = np.array(y_E_C)
y_C_C = np.array(y_C_C)

# 储存数据
np.save(mkpath+'/total_sigma_xx.npy', y_total)
np.save(mkpath+'/B-B_sigma_xx.npy', y_B_B)
np.save(mkpath+'/B-C_sigma_xx.npy', y_B_C)
np.save(mkpath+'/B-E_sigma_xx.npy', y_B_E)
np.save(mkpath+'/E-E_sigma_xx.npy', y_E_E)
np.save(mkpath+'/E-C_sigma_xx.npy', y_E_C)
np.save(mkpath+'/C-C_sigma_xx.npy', y_C_C)

# 同时保存txt格式
np.savetxt(mkpath+'/total_sigma_xx.txt', np.real(y_total))
np.savetxt(mkpath+'/B-B_sigma_xx.txt', np.real(y_B_B))
np.savetxt(mkpath+'/B-C_sigma_xx.txt', np.real(y_B_C))
np.savetxt(mkpath+'/B-E_sigma_xx.txt', np.real(y_B_E))
np.savetxt(mkpath+'/E-E_sigma_xx.txt', np.real(y_E_E))
np.savetxt(mkpath+'/E-C_sigma_xx.txt', np.real(y_E_C))
np.savetxt(mkpath+'/C-C_sigma_xx.txt', np.real(y_C_C))

print("计算完成！")
print(f"总采样点数: {len(omega_range)}")
print(f"保存路径: {mkpath}")
