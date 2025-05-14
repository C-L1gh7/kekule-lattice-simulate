#Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
from math import pi

# t1:inner hopping, t2:inter hppping
t1 = 1
t2 = 2
POE=0.0
NOE=0.0

main_L = 1.0 # length of main lattice
bond_L = 1.0 # bond length between lattices
modulus = 2*main_L + bond_L

edgelength = 25.5

# create a lattice in format of pyBinding
def O_keku(onsite_energy=[POE,POE,POE,NOE,NOE,NOE]):
    "创建一个O型kekule晶格"
    lat = pb.Lattice(
        a1 = [modulus, 0.00000000, 0.00000000],
        a2 = [-modulus/2, np.sqrt(3)/2*modulus, 0.00000000],
        # a3 = [0.00000000, 0.00000000, 10.00000000]
    )
    o = np.array([1/2*main_L, np.sqrt(3)/2*main_L, 0.00000000])
    lat.add_sublattices(
        #name and positions
        ('C1', o + [0.00000000, 0.00000000, 0.00000000], onsite_energy[0]),
        ('C2', o + [1/2*main_L+bond_L, -np.sqrt(3)/2*main_L, 0.00000000], onsite_energy[1]),
        ('C3', o + [-(main_L+1/2*bond_L), np.sqrt(3)/2*bond_L, 0.00000000], onsite_energy[2]),
        ('C4', o + [1/2*main_L, -np.sqrt(3)/2*main_L, 0.00000000], onsite_energy[3]),
        ('C5', o + [1/2*bond_L, np.sqrt(3)/2*bond_L, 0.00000000], onsite_energy[4]),
        ('C6', o + [-main_L, 0.00000000, 0.00000000], onsite_energy[5]),
    )

    lat.add_hoppings(

        #between main cell and the cell (1,1,0)

        #between main cell and the cell (1,0,0)
        ([1, 0, 0], 'C2', 'C6', t1),
        ([1, 0, 0], 'C5', 'C3', t1),

        #between main cell and the cell (1,-1,0)

        #between main cell and the cell (0,1,0)
        ([0, 1, 0], 'C3', 'C4', t1),
        ([0, 1, 0], 'C5', 'C2', t1),

        #inside the main cell
        ([0, 0, 0], 'C1', 'C4', t1),
        ([0, 0, 0], 'C1', 'C5', t2),
        ([0, 0, 0], 'C1', 'C6', t1),
        ([0, 0, 0], 'C2', 'C4', t2),
        ([0, 0, 0], 'C3', 'C6', t2)
    )
    return lat

def trapezoid(o=[0,0], r=10):
    "裁切一个正六边形"
    p0 = np.array(o)# 转换为numpy数组
    p1 = p0 + [1/2*r,np.sqrt(3)/2*r]
    p2 = p0 + [r,0]
    p3 = p0 + [1/2*r,-np.sqrt(3)/2*r]
    p4 = p0 + [-1/2*r,-np.sqrt(3)/2*r]
    p5 = p0 + [-r,0]
    p6 = p0 + [-1/2*r,np.sqrt(3)/2*r]
    return pb.Polygon([p1, p2, p3, p4, p5, p6])

def mkdir(path):
    # 引入模块
    import os

    # 去除首位空格
    path=path.strip()

    # 去除尾部 \ 符号
    path=path.rstrip("\\")

     # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists=os.path.exists(path)

    # 判断结果
    if not isExists:
        # 如果不存在则创建目录
        # 创建目录操作函数
        os.makedirs(path) 
        print (path+' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print (path+' 目录已存在')
        return False 
    
def probability(model,mkpath):
    "计算波函数"
    solver = pb.solver.arpack(model,k=model.system.num_sites-1)
    eigenvalues = solver.calc_eigenvalues()  # position in [nm]
    eigenvalues.plot(mark_degenerate = False, s=0.8,show_indices=True)
    # plt.show()
    plt.savefig(mkpath+'/energy_band_1.eps')  # 保存为PNG文件
    plt.clf()
    solver = pb.solver.arpack(model,k=20)
    eigenvalues = solver.calc_eigenvalues()  # position in [nm]
    eigenvalues.plot(show_indices=True)
    plt.savefig(mkpath+'/energy_band_2.eps')  # 保存为PNG文件
    plt.clf()
    # for i in range(2495,2600):
    #     probability_map = solver.calc_probability([i])
    #     probability_map.plot(site_radius=(0.0, 0.3),cmap=['red'])
    #     plt.savefig(mkpath+'/wave_function/'+str(i)+'.png', dpi=300)  # 保存为PNG文件
    #     plt.clf()
    print("wawe function")


###目录创建函数
def mkdir(path):
    # 引入模块
    import os

    # 去除首位空格
    path=path.strip()

    # 去除尾部 \ 符号
    path=path.rstrip("\\")

     # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists=os.path.exists(path)

    # 判断结果
    if not isExists:
        # 如果不存在则创建目录
        # 创建目录操作函数
        os.makedirs(path) 
        print (path+' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print (path+' 目录已存在')
        return False 