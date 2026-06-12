import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle, FancyArrowPatch

import basic_function as bf

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

####################################################################################
# 全局参数
t1 = 1
t2 = 2 * t1

main_L = 1.0
bond_L = 1.0
modulus = 2 * main_L + bond_L

####################################################################################
# 分别控制三个子图的 L (edgelength) 和 Mu (化学势)
# 子图 (a)
L_a = 49.1       # edgelength = 4.1+3n (有corner态)
Mu_a = 0

# 子图 (b)
L_b = 49.1       # edgelength = 4.1+3n (有corner态)
Mu_b = 0.47

# 子图 (c)
L_c = 48.1       # edgelength ≠ 4.1+3n (无corner态)
Mu_c = 0

####################################################################################
# 数据文件路径（根据 L 值对应不同的数据文件）
# 如果需要重新计算，将 RECOMPUTE 设为 True
RECOMPUTE = False

data_files = {
    48.1: 'data/matrix_output1.txt',
    49.1: 'data/matrix_output2.txt',
}

####################################################################################
# color
cornercolor = '#C71F2D'
bulkcolor = "#38557E"
edgecolor = '#DBA972'
mu_color = '#2CA02C'

s = 0.7

# 控制图的尺寸
fig_width = 8.3        # 总图宽度 (inches)
fig_height = 5        # 总图高度 (inches)
width_ratios = [1, 1, 1]  # 三个子图的宽度比例

mkpath = "result/picture"
bf.mkdir(mkpath)

####################################################################################

def compute_or_load_eigenvalues(L, data_file, recompute=False):
    """计算或加载特征值"""
    lattice = bf.O_keku(t1=t1, t2=t2)
    model = pb.Model(lattice, bf.trapezoid(r=L))

    if recompute or not os.path.exists(data_file):
        solver = pb.solver.arpack(model, k=model.system.num_sites - 1)
        eigenvalues = solver.eigenvalues
        np.savetxt(data_file, eigenvalues)
        print(f"计算完成: L={L}, 保存至 {data_file}")
    else:
        eigenvalues = np.genfromtxt(data_file, dtype=None)
        print(f"加载完成: L={L}, 来自 {data_file}")

    return model, eigenvalues


def classify_states(eigenvalues):
    """将态分类为 corner, edge, bulk"""
    corner_state = int((len(eigenvalues) / 2) - 2)
    values_corner = eigenvalues[corner_state:corner_state + 6]
    CornerState_index = list(range(corner_state, corner_state + 6))

    EdgeState_index = np.where(np.abs(eigenvalues) <= 1)[0].tolist()
    BulkState_index = np.where(np.abs(eigenvalues) > 1)[0].tolist()
    EdgeState = eigenvalues[EdgeState_index]
    BulkState = eigenvalues[BulkState_index]

    return {
        'corner_index': CornerState_index,
        'corner_values': values_corner,
        'edge_index': EdgeState_index,
        'edge_values': EdgeState,
        'bulk_index': BulkState_index,
        'bulk_values': BulkState,
        'corner_state_start': corner_state,
    }


def has_corner_states(L):
    """判断该 edgelength 是否有 corner 态 (edgelength = 4.1 + 3n)"""
    return abs((L - 4.1) % 3) < 0.01


def plot_panel(ax, fig, model, eigenvalues, L, Mu, label, show_ylabel=True, show_inset=True, show_legend=True):
    """绘制单个子图"""
    states = classify_states(eigenvalues)
    num_sites = model.system.num_sites

    if has_corner_states(L):
        ax.scatter(states['edge_index'], states['edge_values'], s=s, color=edgecolor, label='Edge')
        ax.scatter(states['corner_index'], states['corner_values'], s=s, color=cornercolor, label='Corner')
        ax.scatter(states['bulk_index'], states['bulk_values'], s=s, color=bulkcolor, label='Bulk')
    else:
        states_num = list(range(len(eigenvalues)))
        ax.scatter(states_num, eigenvalues, s=0.1, color=bulkcolor)

    ax.set_ylim(-2.5, 2.5)
    ax.set_xlim((num_sites / 2) - 1700, (num_sites / 2) + 1700)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax.tick_params(axis='y', labelsize=12, direction='in')
    ax.set_xlabel('Energy Level', fontsize=14)

    if show_ylabel:
        ax.set_ylabel(r'$\mathrm{Energy}/t_1$', fontsize=14)
    else:
        ax.tick_params(axis='y', which='both', labelleft=False)

    ax.text(-0.1, 1.05, label, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')

    # 化学势虚线
    ax.axhline(y=Mu, color=mu_color, linestyle='--', linewidth=1.0, zorder=0)
    ax.annotate(rf'$\mu={Mu}t_1$', xy=(ax.get_xlim()[1], Mu),
                xytext=(-5, 3), textcoords='offset points',
                fontsize=11, color=mu_color, va='bottom', ha='right', fontweight='bold')

    # inset 图（仅对有 corner 态的子图）
    if has_corner_states(L) and show_inset:
        corner_state = states['corner_state_start']
        axins = ax.inset_axes((0.65, 0.2, 0.3, 0.15))
        axins.scatter(states['edge_index'], states['edge_values'], s=s, color=edgecolor)
        axins.scatter(states['corner_index'], states['corner_values'], s=s, color=cornercolor)
        axins.scatter(states['bulk_index'], states['bulk_values'], s=s, color=bulkcolor)
        axins.set_xlim(corner_state - 1, corner_state + 6)
        axins.set_ylim(-0.0015, 0.0015)
        axins.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        axins.tick_params(axis='y', labelsize=8, direction='in')
        axins.yaxis.set_major_locator(plt.MaxNLocator(3))

        x1, x2 = corner_state - 100, corner_state + 100
        y1, y2 = -5e-2, 5e-2
        inset_indicator_color = '#555555'
        rect = Rectangle((x1, y1), x2 - x1, y2 - y1, edgecolor=inset_indicator_color,
                         facecolor='none', linewidth=0.7, linestyle='--')
        ax.add_patch(rect)

        # 从虚线框到 inset 的箭头
        fig.canvas.draw()
        p1_data = (x2, y1)
        p1_display = ax.transData.transform(p1_data)
        p1_axes = ax.transAxes.inverted().transform(p1_display)
        p2_axes = (0.65, 0.2 + 0.15)  # inset 左上角

        dx = p2_axes[0] - p1_axes[0]
        dy = p2_axes[1] - p1_axes[1]
        length = np.sqrt(dx**2 + dy**2)
        dx_norm = dx / length
        dy_norm = dy / length

        gap_start = 0.03
        gap_end = 0.02
        arrow_start = (p1_axes[0] + gap_start * dx_norm, p1_axes[1] + gap_start * dy_norm)
        arrow_end = (p2_axes[0] - gap_end * dx_norm, p2_axes[1] - gap_end * dy_norm)

        arrow = FancyArrowPatch(
            arrow_start, arrow_end,
            transform=ax.transAxes,
            arrowstyle='->,head_length=3,head_width=2',
            color=inset_indicator_color, linewidth=0.9, linestyle='-'
        )
        ax.add_patch(arrow)

    # 图例
    if has_corner_states(L) and show_legend:
        handles, labels_leg = ax.get_legend_handles_labels()
        new_order = [1, 0, 2]
        handles = [handles[i] for i in new_order]
        labels_leg = [labels_leg[i] for i in new_order]
        ax.legend(handles, labels_leg, frameon=False, markerscale=3)


####################################################################################
# 主程序

# 加载/计算数据
panels = [
    (L_a, Mu_a, '(a)'),
    (L_b, Mu_b, '(b)'),
    (L_c, Mu_c, '(c)'),
]

fig, axs = plt.subplots(1, 3, figsize=(fig_width, fig_height),
                        gridspec_kw={'width_ratios': width_ratios})

show_legend_list = [True, False, False]  # 只有图(a)显示图例

for i, (L, Mu, label) in enumerate(panels):
    data_file = data_files.get(L, f'matrix_output_L{L}.txt')
    model, eigenvalues = compute_or_load_eigenvalues(L, data_file, recompute=RECOMPUTE)
    plot_panel(axs[i], fig, model, eigenvalues, L, Mu, label,
              show_ylabel=(i == 0),
              show_inset=True,
              show_legend=show_legend_list[i])

plt.tight_layout()
output_pdf = os.path.abspath(os.path.join(mkpath, 'energy_band_3panels.pdf'))
plt.savefig(output_pdf)
print(f"保存至: {output_pdf}")

try:
    os.startfile(output_pdf)
except OSError as e:
    print(f"无法自动打开文件: {e}")
