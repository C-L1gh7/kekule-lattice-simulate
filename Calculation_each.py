
# Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from tqdm import tqdm
from numba import njit

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
Mu = 1e-5
Gamma = 0.03
h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength ** 2.
hS = -h_bar / S

####################################################################################

#设置存储位置
mkpath="result/anti-corner"+str(edgelength)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
bf.mkdir(mkpath)

####################################################################################
lattice = bf.O_keku()
lattice.plot()
plt.show()

model = pb.Model(
    lattice,
    bf.trapezoid(r=edgelength)
)

model.plot()
plt.savefig(mkpath+'/model.eps')  # 保存为PNG文件
plt.show() #检查切割形状是否正确
plt.clf()
bf.probability(model,mkpath)

# 计算速度算符
# np.set_printoptions(threshold=sys.maxsize)
length_x = len(model.system.x)
r_x = np.zeros((length_x, length_x), dtype=float)
np.fill_diagonal(r_x, model.system.x)
H = model.hamiltonian.todense()
V_x = (r_x @ H) - (H @ r_x)
np.savetxt(mkpath+'/Vx.txt', V_x)

length_y = len(model.system.y)
r_y = np.zeros((length_y, length_y), dtype=float)
np.fill_diagonal(r_y, model.system.y)
V_y = (r_y @ H) - (H @ r_y)
np.savetxt(mkpath+'/Vy.txt', V_y)

solver = pb.solver.lapack(model)
values = solver.eigenvalues  # 求eigenvalues
vectors = solver.eigenvectors.astype(np.float64)
####################################################################################

indices_to_keep = np.where(np.abs(values) <= 2.5) # 计算的能量范围
values = values[indices_to_keep]
vectors = vectors[indices_to_keep]


corner_state = int((len(values)/2) - 3)
values_corner = values[corner_state:corner_state+6]
CornerState_index = list(range(corner_state, corner_state+6))
plt.scatter(CornerState_index,values[CornerState_index])
# plt.show()
EdgeState_index = np.where(np.abs(values) <= 1)
BulkState_index = np.where(np.abs(values) > 1)
EdgeState_index = EdgeState_index[0].tolist()
BulkState_index = BulkState_index[0].tolist()
plt.scatter(EdgeState_index,values[EdgeState_index])
# plt.show()
plt.scatter(BulkState_index,values[BulkState_index])
# plt.show()
corner_state = int((model.system.num_sites/2) - 3)
values_corner = values[corner_state:corner_state+6]
cornerState_index = list(range(corner_state, corner_state+6))

####################################################################################
@njit
def conductivity(h_omega, CornerState_index, EdgeState_index, BulkState_index): # 计算总光电导以及各态贡献的光电导时使用
# def conductivity(h_omega): # 只计算总光电导的时候使用
    """calculate sigma"""
    # 初始化数值
    sigma = 0
    sigma_B_B = 0
    sigma_C_C = 0
    sigma_B_C = 0
    sigma_B_E = 0
    sigma_E_E = 0
    sigma_E_C = 0
    for n in range(len(vectors)):  # eigenvectors有多长就取多长
        fF_n = 1 / (np.exp((values[n] - Mu) / kBT) + 1)
        aA = np.dot(vectors[n].T, V_x)  # 计算aA = <a|A
        Aa = np.dot(V_x, vectors[n])  # 计算Aa = |A|a>
        for k in range(len(vectors)):
            if n == k:  # bra ket不取同一个eigenvector进行计算
                continue
            else:
                # |a> = φn, |b> = φk, A = V_x
                inner_product_1 = np.dot(aA, vectors[k])  # 计算 <a|A|b>
                inner_product_2 = np.dot(vectors[k].T, Aa)  # 对<b|进行复共轭之后计算 <b|A|a>
                inner_product = inner_product_1 * inner_product_2  # <a|A|b> * <b|A|a>
                fF = fF_n - 1 / (np.exp((values[k] - Mu) / kBT) + 1)  # fF = fF_n - fF_k
                epsilon = values[n] - values[k]
                result_sigma = complex((fF * inner_product.item()), 0) / complex(epsilon * (h_omega + epsilon), Gamma)
                sigma += result_sigma # total
                if n in BulkState_index and k in BulkState_index: # bulk-bulk(B-B)
                    sigma_B_B += result_sigma
                elif (n in BulkState_index and k in CornerState_index) or (n in CornerState_index and k in BulkState_index): # bulk-corner(B-C)
                    sigma_B_C += result_sigma
                elif (n in BulkState_index and k in EdgeState_index and k not in CornerState_index) or (n in EdgeState_index and n not in CornerState_index and k in BulkState_index): # bulk-edge(B-E)
                    sigma_B_E += result_sigma
                elif n in EdgeState_index and k in EdgeState_index and n not in CornerState_index and k not in CornerState_index : # edge-edge(E-E)
                    sigma_E_E += result_sigma
                elif (n in EdgeState_index and k in CornerState_index and n not in CornerState_index) or (n in CornerState_index and k in EdgeState_index and k not in CornerState_index ): # edge-corner(E-C)
                    sigma_E_C +=result_sigma
                elif n in CornerState_index and k in CornerState_index: # corner-corner(C-C)
                    sigma_C_C += result_sigma
                   
    sigma *= complex(0, hS)
    sigma_B_B *= complex(0,hS)
    sigma_B_C *= complex(0,hS)
    sigma_B_E *= complex(0,hS)
    sigma_E_E *= complex(0,hS)
    sigma_E_C *= complex(0,hS)
    sigma_C_C *= complex(0,hS)

    return sigma, sigma_B_B, sigma_B_C, sigma_B_E, sigma_E_E, sigma_E_C, sigma_C_C
    # return sigma   # only total


# 创建一个新的图形
plt.figure()

# 分别绘制全部态和角态参与的光电导
x = list(i for i in np.arange(0,2.505,0.005))
np.savetxt(mkpath+'/h_omega.txt', x)
y_total=[]
# only total
# for i in tqdm(x):
#     conductivity_value_1 = conductivity(i)
#     y_total.append(conductivity_value_1)

y_B_B=[]
y_B_C = []
y_B_E = []
y_E_E = []
y_E_C = []
y_C_C = []
for i in tqdm(x):
    conductivity_value_1, conductivity_value_2, conductivity_value_3, conductivity_value_4, conductivity_value_5, conductivity_value_6, conductivity_value_7 = conductivity(i,CornerState_index, EdgeState_index, BulkState_index)
    y_total.append(conductivity_value_1)
    y_B_B.append(conductivity_value_2)
    y_B_C.append(conductivity_value_3)
    y_B_E.append(conductivity_value_4)
    y_E_E.append(conductivity_value_5)
    y_E_C.append(conductivity_value_6)
    y_C_C.append(conductivity_value_7)

# 储存数据
np.savetxt(mkpath+'/total_sigma_xx.txt', y_total)
np.savetxt(mkpath+'/B-B_sigma_xx.txt', y_B_B)
np.savetxt(mkpath+'/B-C_sigma_xx.txt', y_B_C)
np.savetxt(mkpath+'/B-E_sigma_xx.txt', y_B_E)
np.savetxt(mkpath+'/E-E_sigma_xx.txt', y_E_E)
np.savetxt(mkpath+'/E-C_sigma_xx.txt', y_E_C)
np.savetxt(mkpath+'/C-C_sigma_xx.txt', y_C_C)

# data = [{'y': y_total},{'y': y_B_B},{'y': y_B_C},{'y': y_B_E},{'y': y_E_E},{'y': y_E_C},{'y': y_C_C}]
# # y_max = max(abs(data[i]['y']) for i in range(len(data))) # 设置最大y轴范围
# fig, axs = plt.subplots(4, 2,sharex=True, sharey=True)# 合并x，y轴
# axs = axs.ravel()
# axs[0].plot(x, data[0]['y'], label='B-B', color=(55/255,103/255,149/255), linestyle='-', lw=0.8)
# axs[1].plot(x, data[0]['y'], label='B-C', color=(114/255,188/255,213/255), linestyle='--', lw=0.8)
# axs[2].plot(x, data[0]['y'], label='B-E', color=(255/255,208/255,111/255), linestyle='-.', lw=0.8)
# axs[3].plot(x, data[0]['y'], label='E-E', color=(231/255,98/255,84/255), linestyle=':', lw=0.8)
# axs[4].plot(x, data[0]['y'], label='E-C', color=(55/255,103/255,149/255), linestyle=':', lw=0.8)
# axs[5].plot(x, data[0]['y'], label='B-C', color=(114/255,188/255,213/255), linestyle='--', lw=0.8)
# for i, ax in enumerate(axs):
#     ax.set_aspect(aspect=50.0)
#     ax.set_xlim(0, 2)  # 设置x轴范围
#     # ax.set_ylim(0, y_max)    # 设置y轴范围
# plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05, hspace=0.0)
# plt.savefig(mkpath+'/Gamma_sigma_xx_all.eps')  # 保存为eps文件
# # 绘制曲线
# plt.plot(x, y_total, color= (55/255,103/255,149/255), linestyle='-', lw=0.8, label = 'Total')
# plt.plot(x, y_B_B, color = (114/255,188/255,213/255), linestyle='--', lw=0.8, label = 'B-B')
# plt.plot(x, y_C_C, color = (255/255,208/255,111/255), linestyle='-.', lw=0.8, label = 'C-C')
# plt.plot(x, y_B_C, color = (231/255,98/255,84/255), linestyle=':', lw=0.8, label = 'B-C')
# # 添加标题和轴标签
# plt.legend(frameon=False) # 显示图例
# plt.xlabel(r'$\hbar\omega\ [t]$', fontsize=14)
# plt.ylabel(r'$\sigma_{xx}\ [{\tilde t}^2 e^2/\hbar]$', fontsize=14)
# plt.savefig(mkpath+'/sigma_xx.eps')  # 保存为eps文件