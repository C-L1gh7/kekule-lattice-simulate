# Import libraries
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import math
from tqdm import tqdm
from matplotlib.ticker import AutoMinorLocator
import itertools
import winsound

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

# ==================== 参数设置区域 ====================
# 定义参数数组 - 可以修改这里的值来批量计算不同的参数组合
Mu = np.array([0.2])           # 化学势数组
kBT = np.array([0.01])          # 温度数组
Gamma = np.array([0.001])       # 展宽参数数组
# ====================================================

h_bar = 1.0
S = (3 * math.sqrt(3) / 2) * edgelength ** 2.
hS = -h_bar / S

####################################################################################

# 创建分层存储目录
edge_path = f"result/anti-corner{edgelength}"
bf.mkdir(edge_path)

# 生成所有参数组合
param_combinations = list(itertools.product(Mu, kBT, Gamma))
total_combinations = len(param_combinations)

print(f"总共需要计算 {total_combinations} 组参数组合")

# 设置全局字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

# 降采样函数
def decimate_data(x, y, max_points=5000):
    if len(x) <= max_points:
        return x, y
    step = max(len(x) // max_points, 1)
    indices = np.arange(0, len(x), step)
    return x[indices], y[indices]

# 遍历所有参数组合
for param_idx, (mu_val, kbt_val, gamma_val) in enumerate(param_combinations):
    print(f"\n{'='*60}")
    print(f"正在计算第 {param_idx+1}/{total_combinations} 组: Mu={mu_val}, kBT={kbt_val}, Gamma={gamma_val}")
    print(f"{'='*60}")

    # 设置当前参数组合的存储位置
    mkpath = f"{edge_path}/mu={mu_val}kBT={kbt_val}gamma={gamma_val}"
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
    fF = 1.0 / (np.exp((values - mu_val) / kbt_val) + 1.0)

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
        denominator = epsilon * (h_omega + epsilon) + 1j * gamma_val

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

    for omega in tqdm(omega_range, desc=f"计算各态光电导 (Mu={mu_val}, kBT={kbt_val}, Gamma={gamma_val})"):
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

    ####################################################################################
    # 绘图（参考 all conductivity.py 的风格）
    ####################################################################################

    # 颜色配置（与 all conductivity.py 一致）
    total_color = "#FF6CF3"
    total_color_Between = "#FFB7F9"

    bb_Between = "#8AC7FF"; cb_Between = "#F2DBFD"; eb_Between = "#92FFF0"
    ee_Between = "#9DFFBA"; ce_Between = "#FFFD72"
    bb = "#40A3FF"; cb = "#CC5FFF"; eb = "#24D8C0"
    ee = "#1BA846"; ce = "#D6C800"

    alpha = 0.3
    lw = 0.3
    legend_fontsize = 10
    legend_lw = 1

    # 取实部用于绘图
    x_plot = omega_range
    y_total_real = np.real(y_total)
    y_BB_real = np.real(y_B_B)
    y_BC_real = np.real(y_B_C)
    y_BE_real = np.real(y_B_E)
    y_EE_real = np.real(y_E_E)
    y_EC_real = np.real(y_E_C)

    # 降采样
    x_t,  y_t  = decimate_data(x_plot, y_total_real)
    x_bb_, y_bb_ = decimate_data(x_plot, y_BB_real)
    x_bc_, y_bc_ = decimate_data(x_plot, y_BC_real)
    x_be_, y_be_ = decimate_data(x_plot, y_BE_real)
    x_ee_, y_ee_ = decimate_data(x_plot, y_EE_real)
    x_ec_, y_ec_ = decimate_data(x_plot, y_EC_real)

    # 三个子图
    fig, ax = plt.subplots(3, 1, sharex=True)

    # ===== 图 (a)：Total sigma_xx =====
    ax[0].plot(x_t, y_t, linestyle='-', lw=lw, color=total_color, label='Total')
    ax[0].fill_between(x_t, y_t, color=total_color_Between, alpha=alpha)
    leg0 = ax[0].legend(handlelength=2, frameon=False, loc='upper left', fontsize=legend_fontsize)
    for legline in leg0.get_lines():
        legline.set_linewidth(legend_lw)
    ax[0].tick_params(axis='both', which='both', top=True, labelbottom=False, right=True, direction='in', width=0.5)
    ax[0].set_xlim(0, 2.2)

    # ===== 图 (b)：C-B, C-E, E-E =====
    ax[1].plot(x_bc_, y_bc_, linestyle=':', lw=lw, color=cb, label='C-B')
    ax[1].fill_between(x_bc_, y_bc_, color=cb_Between, alpha=alpha)
    ax[1].plot(x_ec_, y_ec_, linestyle=':', lw=lw, color=ce, label='C-E')
    ax[1].fill_between(x_ec_, y_ec_, color=ce_Between, alpha=alpha)
    ax[1].plot(x_ee_, y_ee_, linestyle='-', lw=lw, color=ee, label='E-E')
    ax[1].fill_between(x_ee_, y_ee_, color=ee_Between, alpha=alpha)
    leg1 = ax[1].legend(fontsize=legend_fontsize, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
    for legline in leg1.get_lines():
        legline.set_linewidth(legend_lw)
    ax[1].minorticks_on()
    ax[1].set_xlim(0, 2.2)
    ax[1].tick_params(axis='both', which='both', top=True, labelbottom=False, right=True, direction='in', width=0.5)
    ax[1].set_ylabel(r'$\mathrm{Re}(\sigma_{xx})/(\tilde{t}^2 e^2/\hbar)$', fontsize=14)

    # ===== 图 (c)：B-B, E-B =====
    ax[2].plot(x_bb_, y_bb_, linestyle='--', lw=lw, color=bb, label='B-B')
    ax[2].fill_between(x_bb_, y_bb_, color=bb_Between, alpha=alpha)
    ax[2].plot(x_be_, y_be_, linestyle=':', lw=lw, color=eb, label='E-B')
    ax[2].fill_between(x_be_, y_be_, color=eb_Between, alpha=alpha)
    leg2 = ax[2].legend(fontsize=legend_fontsize, labelspacing=0.4, handlelength=2, frameon=False, loc='upper left')
    for legline in leg2.get_lines():
        legline.set_linewidth(legend_lw)
    ax[2].minorticks_on()
    ax[2].set_xlim(0, 2.2)
    ax[2].set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.2])
    ax[2].tick_params(axis='both', which='both', top=True, labelbottom=True, right=True, direction='in', width=0.5)
    ax[2].set_xlabel(r'$\hbar\omega/t_1$', fontsize=14)
    ax[2].xaxis.set_tick_params(labelsize=12)

    # 设置y轴刻度字体大小
    ax[0].yaxis.set_tick_params(labelsize=12)
    ax[1].yaxis.set_tick_params(labelsize=12)
    ax[2].yaxis.set_tick_params(labelsize=12)

    # 统一刻度样式
    for a in ax:
        a.tick_params(axis='both', which='major', length=3, width=0.4, direction='in')
        a.tick_params(axis='both', which='minor', length=2, width=0.2, direction='in')
        a.xaxis.set_minor_locator(AutoMinorLocator(5))
        a.yaxis.set_minor_locator(AutoMinorLocator(5))

    # 标注 (a)(b)(c)
    ax[0].text(-0.065, 0.9, '(a)', transform=ax[0].transAxes, fontsize=14, fontweight='bold')
    ax[1].text(-0.065, 0.87, '(b)', transform=ax[1].transAxes, fontsize=14, fontweight='bold')
    ax[2].text(-0.065, 0.87, '(c)', transform=ax[2].transAxes, fontsize=14, fontweight='bold')

    # 调整排版
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

    # 输出 PDF 到结果文件夹
    output_pdf = mkpath + '/optical_conductivity.pdf'
    plt.savefig(output_pdf)
    plt.close()
    print(f"绘图已保存: {output_pdf}")

print(f"\n{'='*60}")
print(f"全部计算完成！共完成 {total_combinations} 组参数")
print(f"{'='*60}")

# 程序结束提示音
winsound.Beep(1000, 500)  # 频率1000Hz，持续500ms
