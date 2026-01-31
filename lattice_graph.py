import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
import numpy as np
import os

# Create directory if it doesn't exist
output_dir = 'result/picture'
os.makedirs(output_dir, exist_ok=True)

main_L = 1.0 # length of main lattice
bond_L = 1.0 # bond length between lattices
modulus = 2*main_L + bond_L
s = 150  # increased circle size
lw = 5  # increased line width
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
edge_length1 = 5.0001  # largest (purple)
edge_length2 = 4.0001  # medium (green)
edge_length3 = 3.0001  # smallest (purple)


# vector
a1 = np.array([modulus, 0])
a2 = np.array([-modulus/2, np.sqrt(3)/2*modulus])

# Create three subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Define three hexagon sizes
hexagon_sizes = [edge_length1, edge_length2, edge_length3]  # [5.5, 4.5, 3.5]
hexagon_colors = [hex_color2, hex_color1, hex_color2]  # [purple, green, purple]

left = -width/2
top = -height/2

for idx, (ax, hex_size, hex_color) in enumerate(zip(axes, hexagon_sizes, hexagon_colors)):
    # Create hexagon path for this subplot
    angle = np.pi / 3  # 60度
    vertices = [
        (hex_size * np.cos(i * angle), hex_size * np.sin(i * angle))
        for i in range(6)
    ]
    hexagon_path = mpltPath.Path(vertices)

    # Calculate the range of lattice cells needed based on hexagon size
    # Use the same zoom level for all subplots
    zoom_factor = 2.5
    view_size = hexagon_sizes[0] * zoom_factor / 2

    # Calculate how many cells we need in each direction
    # Each primitive cell has size 'modulus' in a1 direction and sqrt(3)/2*modulus in a2 direction
    cells_needed_x = int(np.ceil(view_size / modulus)) + 2  # +2 for safety margin
    cells_needed_y = int(np.ceil(view_size / (np.sqrt(3)/2*modulus))) + 2

    # Plot lattice only in the needed range
    for i in np.arange(-cells_needed_x, cells_needed_x + 1):
        for j in np.arange(-cells_needed_y, cells_needed_y + 1):
            # Primitive cell
            o = np.array([1/2*main_L, np.sqrt(3)/2*main_L])
            o += i*a1 + j*a2
            C1 = o
            C2 = o + [1/2*main_L+bond_L, -np.sqrt(3)/2*main_L]
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
                a, b = plotting_order[k][2]
                arr = a*a1 + b*a2
                x2, y2 = plotting_order[k][1] + arr
                is_inside1 = hexagon_path.contains_point((x1, y1))
                is_inside2 = hexagon_path.contains_point((x2, y2))
                if is_inside1 and is_inside2:
                    ax.plot([x1, x2], [y1, y2], color=plotting_order[k][3][0], lw=lw, zorder=0)
                else:
                    ax.plot([x1, x2], [y1, y2], color=plotting_order[k][3][1], lw=lw, zorder=0)

            for k in range(6):
                is_inside = hexagon_path.contains_point((C[k][0], C[k][1]))
                if is_inside:
                    ax.scatter(C[k][0], C[k][1], color=circle_color, s=s, zorder=1)
                else:
                    ax.scatter(C[k][0], C[k][1], color=circle_color_lite, s=s, zorder=1)

    # Draw hexagon boundary
    angles = np.linspace(0, 2*np.pi, 7)
    x_hexagon = hex_size * np.cos(angles)
    y_hexagon = hex_size * np.sin(angles)
    ax.fill(x_hexagon, y_hexagon, color='none', edgecolor=hex_color, lw=3.5, alpha=0.7)

    # Add length annotation on the top edge of hexagon
    # Top edge goes from (hex_size * cos(60°), hex_size * sin(60°)) to (hex_size * cos(120°), hex_size * sin(120°))
    # which is (hex_size/2, hex_size*sqrt(3)/2) to (-hex_size/2, hex_size*sqrt(3)/2)
    top_y = hex_size * np.sqrt(3) / 2
    top_x_left = -hex_size / 2
    top_x_right = hex_size / 2

    # Add vertical markers at ends
    marker_height = 0.3
    ax.plot([top_x_left, top_x_left], [top_y - marker_height, top_y + marker_height],
            color='black', lw=2, zorder=2)
    ax.plot([top_x_right, top_x_right], [top_y - marker_height, top_y + marker_height],
            color='black', lw=2, zorder=2)

    # Draw the top edge with double-sided arrow
    # Adjust arrow positions to end exactly at the inner edge of vertical markers
    ax.annotate('', xy=(top_x_right, top_y), xytext=(top_x_left, top_y),
                arrowprops=dict(arrowstyle='<->', lw=2, color='black', shrinkA=0, shrinkB=0))

    # Add L label above the top edge
    ax.text(0, top_y + 0.4, f'$L = {hex_size:.1f}a_0$',
            fontsize=20,
            verticalalignment='bottom',
            horizontalalignment='center',
            fontfamily='serif',
            fontname='Times New Roman',
            bbox=dict(boxstyle='round,pad=0.3',
                     facecolor='white',
                     alpha=0.8,
                     edgecolor='none'))

    # Add a_0 annotations for the first subplot only (idx == 0)
    if idx == 0:
        # Center hexagon coordinates (when i=0, j=0)
        o_center = np.array([1/2*main_L, np.sqrt(3)/2*main_L])
        C1_center = o_center
        C2_center = o_center + [1/2*main_L+bond_L, -np.sqrt(3)/2*main_L]
        C4_center = o_center + [1/2*main_L, -np.sqrt(3)/2*main_L]
        C6_center = o_center + [-main_L, 0]

        # Annotate t1 bond (C6-C1, horizontal red bond)
        t1_y = C1_center[1]
        t1_x_left = C6_center[0]
        t1_x_right = C1_center[0]

        # Draw double-sided arrow for t1
        ax.annotate('', xy=(t1_x_right, t1_y), xytext=(t1_x_left, t1_y),
                    arrowprops=dict(arrowstyle='<->', lw=2, color='black', shrinkA=0, shrinkB=0))

        # Add vertical markers for t1
        ax.plot([t1_x_left, t1_x_left], [t1_y - 0.15, t1_y + 0.15],
                color='black', lw=2, zorder=3)
        ax.plot([t1_x_right, t1_x_right], [t1_y - 0.15, t1_y + 0.15],
                color='black', lw=2, zorder=3)

        # Add a_0 label for t1 (positioned below)
        ax.text((t1_x_left + t1_x_right) / 2, t1_y - 0.35, '$a_0$',
                fontsize=18,
                verticalalignment='top',
                horizontalalignment='center',
                fontfamily='serif',
                fontname='Times New Roman',
                bbox=dict(boxstyle='round,pad=0.2',
                         facecolor='white',
                         alpha=0.8,
                         edgecolor='none'))

        # Annotate t2 bond (C4-C2, horizontal blue bond)
        t2_y = C2_center[1]
        t2_x_left = C4_center[0]
        t2_x_right = C2_center[0]

        # Draw double-sided arrow for t2
        ax.annotate('', xy=(t2_x_right, t2_y), xytext=(t2_x_left, t2_y),
                    arrowprops=dict(arrowstyle='<->', lw=2, color='black', shrinkA=0, shrinkB=0))

        # Add vertical markers for t2
        ax.plot([t2_x_left, t2_x_left], [t2_y - 0.15, t2_y + 0.15],
                color='black', lw=2, zorder=3)
        ax.plot([t2_x_right, t2_x_right], [t2_y - 0.15, t2_y + 0.15],
                color='black', lw=2, zorder=3)

        # Add a_0 label for t2 (positioned below)
        ax.text((t2_x_left + t2_x_right) / 2, t2_y - 0.35, '$a_0$',
                fontsize=18,
                verticalalignment='top',
                horizontalalignment='center',
                fontfamily='serif',
                fontname='Times New Roman',
                bbox=dict(boxstyle='round,pad=0.2',
                         facecolor='white',
                         alpha=0.8,
                         edgecolor='none'))

    # Set plot properties with zoomed-in view
    # Use the same zoom level for all subplots (based on the first/largest hexagon)
    view_width = view_size * 2
    view_height = view_size * 2

    ax.set_xlim(-view_width/2, view_width/2)
    ax.set_ylim(-view_height/2, view_height/2)
    ax.set_aspect('equal')
    ax.axis('off')

    # Add label (a, b, c) in top-left corner with Times New Roman font
    labels = ['a', 'b', 'c']
    label_text = ax.text(0.025, 0.975, f'({labels[idx]})',
            transform=ax.transAxes,
            fontsize=28,
            verticalalignment='top',
            horizontalalignment='left',
            fontfamily='serif',
            fontname='Times New Roman',
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.1',
                     facecolor='white',
                     alpha=0.7,
                     edgecolor='none'))

plt.tight_layout()
plt.savefig('result/picture/lattice_graph.pdf', bbox_inches='tight')
