import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import basic_function as bf

# 设置全局字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

####################################################################################
# 参数设置
main_L = 1.0
bond_L = 1.0
modulus = 2 * main_L + bond_L
edgelength = 49.1  # 六边形边长

# 设置存储位置
mkpath = "result/picture"
bf.mkdir(mkpath)

####################################################################################
# 扫描参数范围
t1_values = np.linspace(0.1, 3.0, 40)  # 增加到40个点
t2_values = np.linspace(0.1, 3.0, 40)

# 存储相图数据
phase_data = np.zeros((len(t2_values), len(t1_values)))

print("开始计算相图...")

# 定义判断阈值（在循环外定义，以便在图例中使用）
zero_energy_threshold = 0.001  # 零能阈值
gap_threshold = 0.1  # 能隙阈值
min_corner_states = 2  # 最少角态数量

for i, t2 in enumerate(tqdm(t2_values, desc="扫描 t2")):
    for j, t1 in enumerate(t1_values):
        try:
            # 创建晶格
            lattice = bf.O_keku(t1=t1, t2=t2)

            # 创建模型
            model = pb.Model(
                lattice,
                bf.trapezoid(r=edgelength)
            )

            # 求解零能附近的本征值（使用 shift-invert 模式）
            num_sites = model.system.num_sites
            k = min(200, num_sites - 2)  # 增加求解能级数量到200

            # 使用 sigma=0 从零能附近开始求解
            solver = pb.solver.arpack(model, k=k, sigma=0)
            eigenvalues = solver.eigenvalues

            # 排序能量本征值
            sorted_E = np.sort(eigenvalues)

            # 找到零能附近的态
            gap_states_mask = np.abs(sorted_E) < zero_energy_threshold
            gap_states = sorted_E[gap_states_mask]

            # 找到能隙位置（正负能量的分界）
            negative_E = sorted_E[sorted_E < -zero_energy_threshold]
            positive_E = sorted_E[sorted_E > zero_energy_threshold]

            # 判断相的标准：
            # 1. SOTI相：存在能隙 + 能隙中有少量零能态（角态）
            # 2. 普通绝缘体：存在能隙但没有零能态
            # 3. 金属相：没有能隙

            if len(negative_E) > 0 and len(positive_E) > 0:
                # 计算能隙大小
                gap = positive_E[0] - negative_E[-1]

                if gap > gap_threshold:  # 存在明显能隙
                    if len(gap_states) >= min_corner_states:  # 至少2个零能态（六边形应该有6个角态）
                        # SOTI 相：有能隙且能隙中有零能态
                        phase_data[i, j] = 1
                    else:
                        # 普通绝缘体：有能隙但没有零能态
                        phase_data[i, j] = 0
                else:
                    # 金属相：能隙很小或不存在
                    phase_data[i, j] = 2
            else:
                # 无法判断，可能是金属相
                phase_data[i, j] = 2

        except Exception as e:
            print(f"\n计算失败 at t1={t1:.2f}, t2={t2:.2f}: {e}")
            phase_data[i, j] = -1  # 标记为计算失败

print("相图计算完成！")

# 保存数据
np.savetxt(mkpath + '/phase_diagram_data.txt', phase_data)
np.savetxt(mkpath + '/t1_values.txt', t1_values)
np.savetxt(mkpath + '/t2_values.txt', t2_values)

####################################################################################
# 绘制相图
fig, ax = plt.subplots(figsize=(10, 8))

# 定义颜色和标签
# 0: NI (Normal Insulator) - 红色
# 1: SOTI (Second-Order Topological Insulator) - 灰色
# 2: Metal - 蓝色
colors = ['#E74C3C', '#95A5A6', '#3498DB']
phase_names = ['Normal Insulator (NI)', 'SOTI', 'Metal']
phase_descriptions = [
    'Gap exists, no corner states',
    'Gap exists + corner states',
    'No gap (gapless)'
]

from matplotlib.patches import Patch

# 绘制相图
T1, T2 = np.meshgrid(t1_values, t2_values)
im = ax.contourf(T1, T2, phase_data, levels=[-0.5, 0.5, 1.5, 2.5],
                 colors=colors, alpha=0.8)

# 标注你的参数点 (t1=1, t2=2)
ax.plot(1, 2, marker='*', color='#FFD700', markersize=20,
        markeredgecolor='black', markeredgewidth=2,
        label=r'This work: $t_1=1, t_2=2$', zorder=10)

# 设置坐标轴
ax.set_xlabel(r"$t_1$ (intra-lattice hopping)", fontsize=18)
ax.set_ylabel(r"$t_2$ (inter-lattice hopping)", fontsize=18)
ax.tick_params(axis='both', labelsize=14, direction='in', width=1.5, length=6)
ax.set_xlim(t1_values[0], t1_values[-1])
ax.set_ylim(t2_values[0], t2_values[-1])

# 创建详细的图例（包含判断条件）
legend_elements = [
    Patch(facecolor=colors[0], edgecolor='black',
          label=f'NI: Normal Insulator\n(Gap > {gap_threshold}, |E| < {zero_energy_threshold}: 0 states)'),
    Patch(facecolor=colors[1], edgecolor='black',
          label=f'SOTI: 2nd-Order TI\n(Gap > {gap_threshold}, |E| < {zero_energy_threshold}: ≥{min_corner_states} states)'),
    Patch(facecolor=colors[2], edgecolor='black',
          label=f'Metal: Metallic phase\n(Gap ≤ {gap_threshold})'),
    plt.Line2D([0], [0], marker='*', color='#FFD700', markersize=12,
               markeredgecolor='black', markeredgewidth=1.5, linestyle='None',
               label=r'Your parameters ($t_1=1, t_2=2$)')
]

ax.legend(handles=legend_elements, fontsize=12, loc='upper right',
          frameon=True, fancybox=True, shadow=True,
          framealpha=0.95, edgecolor='black', facecolor='white')

# 添加网格
ax.grid(True, alpha=0.3, linestyle=':', linewidth=1, color='gray')

# 添加标题
ax.set_title('Phase Diagram of Kekulé Lattice', fontsize=20, fontweight='bold', pad=15)

plt.tight_layout()
plt.savefig(mkpath + '/phase_diagram.pdf', dpi=300, bbox_inches='tight')
plt.savefig(mkpath + '/phase_diagram.png', dpi=300, bbox_inches='tight')
print(f"\n相图已保存到 {mkpath}/phase_diagram.pdf")
print(f"相图已保存到 {mkpath}/phase_diagram.png")

plt.show()
