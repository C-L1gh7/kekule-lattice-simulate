import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
import numpy as np

main_L = 1.0 # length of main lattice
bond_L = 1.0 # bond length between lattices
modulus = 2*main_L + bond_L
s = 20
lw = 1
t1_color = '#CD442F' #red
t2_color = '#13679E' # blue
circle_color = '#D07F2C'# orange

t1_color_lite = '#CC8274' #red
t2_color_lite = '#9BC7D6' # blue
circle_color_lite = '#D1A56A'# orange

hex_color1 = '#33D47E' # green
hex_color2 = '#BD4EFD' # purple


number = 8
# rectangle
width = 40.1
height = 30.5

#hexagon
edge_length1 = 4.5
edge_length2 = 3.5
edge_length3 = 5.5


# vector
a1 = np.array([modulus, 0])
a2 = np.array([-modulus/2, np.sqrt(3)/2*modulus])

angle = np.pi / 3  # 60度
vertices = [
    (edge_length1 * np.cos(i * angle), edge_length1 * np.sin(i * angle))
    for i in range(6)
]

# 创建正六边形的路径对象
hexagon_path = mpltPath.Path(vertices)

o = np.array([1/2*main_L, np.sqrt(3)/2*main_L])
left = -width/2
top = -height/2
for i in np.arange(-number+1,number):
    for j in np.arange(-number+1,number):
        # Primitive cell
        o = np.array([1/2*main_L, np.sqrt(3)/2*main_L])
        o += i*a1 + j*a2
        C1= o
        C2 = o +[1/2*main_L+bond_L, -np.sqrt(3)/2*main_L]
        C3 = o + [-(main_L+1/2*bond_L), np.sqrt(3)/2*bond_L]
        C4 = o + [1/2*main_L, -np.sqrt(3)/2*main_L]
        C5 = o + [1/2*bond_L, np.sqrt(3)/2*bond_L]
        C6 = o + [-main_L, 0]
        C = [C1, C2, C3, C4, C5, C6]

        plotting_order = [
            [C1, C4, [0, 0], [t1_color, t1_color_lite]],
            [C1, C5, [0, 0], [t2_color, t2_color_lite]],
            [C1, C6, [0, 0], [t1_color, t1_color_lite]],
            [C2, C4, [0, 0], [t2_color, t2_color_lite]],
            [C3, C6, [0, 0], [t2_color, t2_color_lite]],
            [C2, C6, [1, 0], [t1_color, t1_color_lite]],
            [C5, C3, [1, 0], [t1_color, t1_color_lite]],
            [C3, C4, [0, 1], [t1_color, t1_color_lite]],
            [C5, C2, [0, 1], [t1_color, t1_color_lite]]
        ]
        for k in range(len(plotting_order)):
            x1, y1 = plotting_order[k][0]
            a, b =plotting_order[k][2]
            arr = a*a1 + b*a2
            x2, y2 = plotting_order[k][1]+ arr
            is_inside1 = hexagon_path.contains_point((x1, y1))
            is_inside2 = hexagon_path.contains_point((x2, y2))
            if is_inside1 and is_inside2:
                plt.plot([x1, x2], [y1, y2], color=plotting_order[k][3][0], lw=lw, zorder=0)
            else:
                plt.plot([x1, x2], [y1, y2], color=plotting_order[k][3][1], lw=lw, zorder=0)

        for k in range(6):
            # print(f"C{k+1}: {C[k]}")
            is_inside = hexagon_path.contains_point((C[k][0], C[k][1]))
            if is_inside:
                plt.scatter(C[k][0], C[k][1], color=circle_color, s=s, zorder=1)
            else:
                plt.scatter(C[k][0], C[k][1], color=circle_color_lite, s=s, zorder=1)

# hexagon
angles = np.linspace(0, 2*np.pi, 7)
x_hexagon = edge_length1 * np.cos(angles)
y_hexagon = edge_length1 * np.sin(angles)
plt.fill(x_hexagon, y_hexagon, color='none', edgecolor=hex_color1, lw=2, alpha=0.7)

x_hexagon = edge_length2 * np.cos(angles)
y_hexagon = edge_length2 * np.sin(angles)
plt.fill(x_hexagon, y_hexagon, color='none', edgecolor=hex_color2, lw=1.2, alpha=0.7)

x_hexagon = edge_length3 * np.cos(angles)
y_hexagon = edge_length3 * np.sin(angles)
plt.fill(x_hexagon, y_hexagon, color='none', edgecolor=hex_color2, lw=1.2, alpha=0.7)

plt.xlim(left, left+width)
plt.ylim(top+height, top)

plt.gca().set_aspect('equal')
plt.axis('off')
plt.savefig('result/picture/lattice_graph.pdf')
