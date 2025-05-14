#Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
from math import pi

# t1:inner hopping, t2:inter hppping
t1=1.0
t2=2.0
POE=0.0
NOE=0.0

# create a lattice in format of pyBinding
def GetLattice(onsite_energy=[POE,POE,POE,NOE,NOE,NOE]):
    lat = pb.Lattice(
        a1 = [4.27542162, 0.00000000, 0.00000000],
        a2 = [-2.13771081, 3.70262373, 0.00000000],
        # a3 = [0.00000000, 0.00000000, 10.00000000]
    )

    lat.add_sublattices(
        #name and positions
        ('C1', [0.71257061, 1.23420850, 0.00000000], onsite_energy[0]),
        ('C2', [2.85028040, 0.00000000, 0.00000000], onsite_energy[1]),
        ('C3', [-1.42514020, 2.46841523, 0.00000000], onsite_energy[2]),
        ('C4', [1.42514122, 0.00000000, 0.00000000], onsite_energy[3]),
        ('C5', [1.42514020, 2.46841523, 0.00000000], onsite_energy[4]),
        ('C6', [-0.71257061, 1.23420850, 0.00000000], onsite_energy[5]),
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

# def make_k_path(k1, k2, step=0.01, **kwargs):
#     #either choose num_steps or step
#     num_steps = 1
#     if 'num_steps' in kwargs:
#         num_steps = kwargs['num_steps']
#     else:
#         num_steps = int(np.linalg.norm(k2 - k1) // step)

#     #k_path.shape == num_steps, k_space_dimensions
#     k_path = np.array([np.linspace(s, e, num_steps, endpoint=False)
#                        for s, e in zip(k1, k2)]).T
#     return k_path


#setup lattice with on-site potential terms
lattice = GetLattice()
# plt.figure()
# plt.subplot(121)
# plt.title('Lattice: xy')
# lat.plot()
# plt.savefig('lattice.png')

# plt.subplot(122)
# plt.title('Lattice: yz')
# lat.plot(axes='yz')
# plt.show()





# $$$$$$ band write code:$$$$$$$
# with open('big.dat', 'w') as f:
#     i = 0
#     for k in full_kpath:
#         model.set_wave_vector(k)
#         solver = pb.solver.lapack(model)
#         tmp = solver.eigenvalues
#         f.write(str(i)+' ')
#         f.write(' '.join(map(str, tmp)))
#         f.write('\n')
#         i += 1

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
################


# 截取形状:kekule-armchair
def trapezoid(o=[0,0], r=10):
    p0 = np.array(o)# 转换为numpy数组
    p1 = p0 + [r,0]
    p2 = p0 + [0.5*r, -1*np.sqrt(3)*r/2]
    p3 = p0 + [-0.5*r, -1*np.sqrt(3)*r/2]
    p4 = p0 + [-1*r, 0.0]
    p5 = p0 + [-0.5*r, np.sqrt(3)*r/2]
    p6 = p0 + [0.5*r, np.sqrt(3)*r/2]
    return pb.Polygon([p1, p2, p3 ,p4 ,p5 ,p6])

edgelength=45.5

model = pb.Model(
    lattice,
    trapezoid(r=edgelength)
    # trapezoid(r=35)
    # trapezoid(r=48)
#     # trapezoid(r=62)
)

# model.plot()
# model.shape.plot()
# plt.savefig(str(edgelength)+'-shape-more.png')



# solver = pb.solver.lapack(model)


mkpath="anti-corner"+str(edgelength)# 打印文件名

# 调用函数
# mkdir(mkpath)# 创建文件夹和文件

solver = pb.solver.arpack(model, k=20)
eigenvalues = solver.calc_eigenvalues()
eigenvalues.plot()
eigenvalues.plot_heatmap(show_indices=True)
plt.show()

for i in range(7,13):
    probability_map = solver.calc_probability([i])
    probability_map.plot(site_radius=(0.0, 0.3),cmap=['red'])
    plt.show()

# # # 分开画角态分布
# for i in range(1392,1407):
#     probability_map = solver.calc_probability([i])
#     probability_map.plot(hopping={'width': 0}, site_radius=(0.0, 1.0), cmap=['red'])
#     ### 上面hopping是指画不画model里头的hopping，site_radius的两个值分别是最小dos和最大所对应的圆的半径
#     plt.savefig(mkpath+'/'+str(edgelength)+'-'+str(i)+'-.png')
#     plt.cla()

# print(model.system.num_sites)

# with open(str(edgelength)+'-eig-more.dat', 'w') as f:
#     i = 0
#     tmp = solver.eigenvalues
#     for item in tmp:
#         # print(item)
#         f.write(str(i)+' '+str(item))
#         f.write('\n')
#         i += 1