
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
t2 = 1.8*t1
POE = 0.0
NOE = 0.0

main_L = 1.0  # length of main lattice
bond_L = 1.0  # bond length between lattices
modulus = 2 * main_L + bond_L
edgelength = 49.1 # edgelength = 4.1+3n(n=0,1,2,...)
####################################################################################

kBT = 0.01
Mu = 0.0
Gamma = 0.001
h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength ** 2.
hS = -h_bar / S

####################################################################################

#设置存储位置
mkpath="result/anti-corner"+str(edgelength)+"_"+str(t2/t1)+"/mu="+str(Mu)+"kBT="+str(kBT)+"gamma="+str(Gamma)# 打印文件名
bf.mkdir(mkpath)

####################################################################################
lattice = bf.O_keku(t1=t1, t2=t2)
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
def conductivity(h_omega): # 只计算总光电导的时候使用
    """calculate sigma"""
    # 初始化数值
    sigma = 0
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
                # if n in BulkState_index and k in BulkState_index: # bulk-bulk(B-B)
                #     sigma_B_B += result_sigma
                # elif (n in BulkState_index and k in CornerState_index) or (n in CornerState_index and k in BulkState_index): # bulk-corner(B-C)
                #     sigma_B_C += result_sigma
                # elif (n in BulkState_index and k in EdgeState_index and k not in CornerState_index) or (n in EdgeState_index and n not in CornerState_index and k in BulkState_index): # bulk-edge(B-E)
                #     sigma_B_E += result_sigma
                # elif n in EdgeState_index and k in EdgeState_index and n not in CornerState_index and k not in CornerState_index : # edge-edge(E-E)
                #     sigma_E_E += result_sigma
                # elif (n in EdgeState_index and k in CornerState_index and n not in CornerState_index) or (n in CornerState_index and k in EdgeState_index and k not in CornerState_index ): # edge-corner(E-C)
                #     sigma_E_C +=result_sigma
                # elif n in CornerState_index and k in CornerState_index: # corner-corner(C-C)
                #     sigma_C_C += result_sigma
                   
    sigma *= complex(0, hS)
    return sigma   # only total


# 创建一个新的图形
plt.figure()

# 分别绘制全部态和角态参与的光电导
x = list(i for i in np.arange(0,0.5,0.005))
np.savetxt(mkpath+'/h_omega.txt', x)
y_total=[]
# only total
for i in tqdm(x):
    conductivity_value_1 = conductivity(i)
    y_total.append(conductivity_value_1)


# 储存数据
np.savetxt(mkpath+'/total_sigma_xx.txt', y_total)