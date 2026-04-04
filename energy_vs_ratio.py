import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pybinding as pb
from tqdm import tqdm

import basic_function as bf


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"


main_L = 1.0
bond_L = 1.0
modulus = 2 * main_L + bond_L
edgelength = 49.1

mkpath = "result/picture"
datapath = "data"
bf.mkdir(mkpath)
bf.mkdir(datapath)


t1_fixed = 1.0
ratio_start = 0
ratio_end = 2.5
ratio_step = 0.005
ratio_values = np.arange(ratio_start, ratio_end + ratio_step, ratio_step)
t2_values = ratio_values * t1_fixed


k_eigenvalues = 80
params = {
    "edgelength": edgelength,
    "t1_fixed": t1_fixed,
    "ratio_start": ratio_start,
    "ratio_end": ratio_end,
    "ratio_step": ratio_step,
    "k_eigenvalues": k_eigenvalues,
}

params_file = os.path.join(datapath, "params.json")
data_file = os.path.join(datapath, "energy_bands.npy")
ratio_file = os.path.join(datapath, "ratio_values.npy")


need_recalculate = True
if os.path.exists(params_file) and os.path.exists(data_file) and os.path.exists(ratio_file):
    with open(params_file, "r", encoding="utf-8") as f:
        saved_params = json.load(f)

    if saved_params == params:
        print("Parameters unchanged, loading cached data...")
        energy_bands = np.load(data_file)
        ratio_values = np.load(ratio_file)
        need_recalculate = False
        print("Cached data loaded.")
    else:
        print("Parameters changed, recalculating...")
else:
    print("No cached data found, starting calculation...")


if need_recalculate:
    energy_data = []
    energy_bands = []
    print("Calculating energy bands versus t2/t1...")

    for ratio, t2 in tqdm(
        zip(ratio_values, t2_values),
        total=len(ratio_values),
        desc="Scanning t2/t1",
    ):
        try:
            lattice = bf.O_keku(t1=t1_fixed, t2=t2)
            model = pb.Model(lattice, bf.trapezoid(r=edgelength))

            num_sites = model.system.num_sites
            k = min(k_eigenvalues, num_sites - 2)

            solver = pb.solver.lapack(model)
            all_eigenvalues = solver.eigenvalues

            sorted_eigenvalues = np.sort(all_eigenvalues)
            num_total = len(sorted_eigenvalues)
            center_idx = num_total // 2
            half_k = k // 2

            start_idx = max(0, center_idx - half_k)
            end_idx = min(num_total, center_idx + half_k)
            eigenvalues = sorted_eigenvalues[start_idx:end_idx]

            for energy in eigenvalues:
                energy_data.append([ratio, energy])

            energy_bands.append(eigenvalues)
        except Exception as e:
            print(f"\nCalculation failed at t2/t1={ratio:.3f}: {e}")

    print("Energy band calculation finished.")

    energy_data = np.array(energy_data)
    energy_bands = np.array(energy_bands)

    np.save(data_file, energy_bands)
    np.save(ratio_file, ratio_values)
    with open(params_file, "w", encoding="utf-8") as f:
        json.dump(params, f, indent=2)
    print(f"Data saved to {datapath}/")


energy_data = []
for i, ratio in enumerate(ratio_values):
    for energy in energy_bands[i]:
        energy_data.append([ratio, energy])
energy_data = np.array(energy_data)


np.savetxt(
    os.path.join(mkpath, "energy_vs_ratio_data.txt"),
    energy_data,
    header="Column 1: t2/t1 ratio, Column 2: Energy",
)


fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(energy_data[:, 0], energy_data[:, 1], s=0.5, c="#134170", alpha=0.8)

ax.axvline(
    x=2,
    color="green",
    linestyle="-",
    linewidth=2.5,
    label=r"This work: $t_2/t_1=2$",
    alpha=0.8,
)

ax.set_xlabel(r"$t_2/t_1$", fontsize=22)
ax.set_ylabel(r"Energy/$t_1$", fontsize=22)
ax.tick_params(axis="both", labelsize=20, direction="in", width=1.5, length=6)
ax.set_xlim(ratio_values[0], ratio_values[-1])

yticks = ax.get_yticks()
yticklabels = [f"{y:.1f}" if y != yticks[0] else "" for y in yticks]
ax.set_yticklabels(yticklabels)

ax.legend(fontsize=20, loc="upper left", frameon=False)

plt.tight_layout()
output_pdf = os.path.abspath(os.path.join(mkpath, "energy_vs_ratio.pdf"))
output_png = os.path.abspath(os.path.join(mkpath, "energy_vs_ratio.png"))

plt.savefig(output_pdf, bbox_inches="tight")
plt.savefig(output_png, dpi=300, bbox_inches="tight")
print(f"\nSaved figure to {output_pdf}")
print(f"Saved figure to {output_png}")

try:
    os.startfile(output_pdf)
except OSError as e:
    print(f"Failed to open output file automatically: {e}")

plt.show()
