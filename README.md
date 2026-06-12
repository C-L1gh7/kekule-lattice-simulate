# Kekulé Lattice Simulation

本项目基于 [pybinding](https://pybinding.site/) 框架构建 Kekulé-type 蜂窝晶格，计算其能带结构、波函数分布与光电导（optical conductivity），用于研究高阶拓扑绝缘体中**角态（corner state）**、**边界态（edge state）** 和**体态（bulk state）** 对光电导的贡献。

## 项目结构

```
.
├── basic_function.py                  # 核心模块：晶格定义、几何裁剪、工具函数
├── Calculation_each_optimized.py      # 主计算脚本：光电导分解（角态/边界态/体态贡献）
├── Calculation_total_optimized.py     # 总光电导计算（无状态分类）
├── EnergyBand_3panels.py              # 三种六边形纳米盘能带对比图
├── wavefunction.py                    # 波函数概率分布可视化
├── phase_diagram.py                   # 相图计算（t1-t2 参数扫描）
├── energy_vs_ratio.py                 # 能量本征值 vs t₂/t₁ 比率扫描
├── opticalConductivity_MuCompare.py   # 化学势 μ 对比图（L=49 / L=48）
├── opticalconductivity_GammaCompare.py # 展宽 Γ 对比图（L=49 / L=48）
├── all conductivity.py                # 光电导带峰位标注的详细图
├── lattice_graph.py                   # 晶格几何示意图
├── LICENSE                            # 许可证
├── data/                              # 缓存的计算数据（.npy / .json）
├── result/                            # 所有输出结果（图片、数据文件）
└── .gitignore
```

## 功能概览

### 核心计算

| 脚本 | 功能 | 输出 |
|------|------|------|
| `Calculation_each_optimized.py` | 构建 Kekulé 晶格模型，对角化哈密顿量，将本征态分类为角态/边界态/体态，计算各类态间跃迁对光电导 σₓₓ(ω) 的贡献 | `result/anti-corner{L}/` 下的能带图 PDF、光电导三面板 PDF、`.npy`/`.txt` 数据文件 |
| `Calculation_total_optimized.py` | 同上但不做状态分类，仅计算总光电导 | `result/edge{L}/` 下的光电导数据 |

### 可视化与分析

| 脚本 | 功能 | 输出 |
|------|------|------|
| `EnergyBand_3panels.py` | 三种六边形纳米盘能带对比：(a) L=49.1 μ=0（有角态）、(b) L=49.1 μ=0.47、(c) L=48.1 μ=0（无角态），带 inset 放大角态区域 | `result/picture/energy_band_3panels.pdf` |
| `wavefunction.py` | 绘制角态与边界态波函数的实空间概率分布 | `result/picture/wavefunction.pdf` |
| `all conductivity.py` | 光电导图：(a) 总光电导 L=49 vs L=48、(b) C-B / C-E / E-E、(c) B-B / E-B，带峰位箭头标注 | `result/picture/optical_conductivity2.pdf` |
| `opticalConductivity_MuCompare.py` | 不同化学势 μ 下光电导对比（纳米盘 L=49 / L=48） | `result/picture/optical_conductivity_MuCompare.pdf` |
| `opticalconductivity_GammaCompare.py` | 不同展宽 Γ 下光电导对比（纳米盘 L=49 / L=48） | `result/picture/optical conductivity GammaCompare.pdf` |
| `lattice_graph.py` | 绘制 Kekulé 晶格几何结构示意图（三种六边形尺寸） | `result/picture/lattice_graph.pdf` |
| `phase_diagram.py` | 扫描 t₁-t₂ 参数空间，分类 SOTI / Normal Insulator / Metal 相 | `result/picture/phase_diagram.pdf` |
| `energy_vs_ratio.py` | 扫描 t₂/t₁ 比率，绘制能带演化图（带数据缓存） | `result/picture/energy_vs_ratio.pdf` |

### 基础模块

[`basic_function.py`](basic_function.py) 提供：
- **`O_keku(t1, t2, onsite_energy)`** — 创建 O-type Kekulé 晶格（6 子晶格原胞）
- **`trapezoid(o, r)`** — 定义正六边形裁剪区域
- **`probability(model, mkpath)`** — 求解并绘制能带图
- **`mkdir(path)`** — 递归创建目录

## 物理模型

- **晶格类型**：Kekulé-O 型变形蜂窝晶格，原胞含 6 个子晶格
- **跃迁参数**：t₁（原胞内跃迁）和 t₂（原胞间跃迁），默认 t₁ = 1, t₂ = 2
- **裁剪形状**：正六边形，边长 `edgelength = 4.1 + 3n` 时支持角态（高阶拓扑绝缘体相）
- **光电导公式**：基于 Kubo 公式，计算实部 Re[σₓₓ(ω)]

## 环境依赖

- Python 3.8.5（考虑到 pybinding 对较新版本 Python 的支持限制）
- [pybinding](https://pybinding.site/) 及其 C++ 编译依赖
- `numpy` `matplotlib` `tqdm` `numba`

## 快速开始

1. **运行主计算**：
   ```bash
   python Calculation_each_optimized.py
   ```
   该脚本会构建六边形 Kekulé 晶格、对角化并计算光电导分解。可在脚本顶部的参数区域修改 `edgelength`、`Mu`、`kBT`、`Gamma` 等。

2. **生成论文用图**：
   ```bash
   python EnergyBand_3panels.py          # 能带对比图
   python all conductivity.py             # 光电导详细图
   python opticalConductivity_MuCompare.py # μ 依赖对比
   python opticalconductivity_GammaCompare.py # Γ 依赖对比
   python wavefunction.py                 # 波函数分布图
   python lattice_graph.py                # 晶格示意图
   python phase_diagram.py                # 相图
   python energy_vs_ratio.py              # 能带 vs t₂/t₁
   ```

3. **查看结果**：所有输出保存在 `result/` 目录下，包括 PDF 矢量图和 `.npy`/`.txt` 数据文件。

## 关键参数

| 参数 | 含义 | 默认值 | 位置 |
|------|------|--------|------|
| `t1` | 原胞内跃迁振幅 | 1 | 各脚本顶部 |
| `t2` | 原胞间跃迁振幅 | 2 | 各脚本顶部 |
| `edgelength` | 六边形边长（控制系统尺寸和拓扑相） | 16 / 49 | 各脚本顶部 |
| `Mu` / `μ` | 化学势 | 0 | 参数设置区域 |
| `kBT` | 温度（能量单位） | 0.01 | 参数设置区域 |
| `Gamma` / `Γ` | 展宽参数 | 0.001 | 参数设置区域 |

## 输出说明

- **`result/anti-corner{L}/`** — 含角态系统（L = 4a_0 + 3na_0）的计算结果
- **`result/edge{L}/`** — 无角态系统（L ≠ 4a_0 + 3na_0）的计算结果
- **`result/picture/`** — 汇总对比图和示意图
- 每个子目录包含：能带图、光电导图（PDF）、本征值/本征矢/速度算符/光电导数据（`.npy` + `.txt`）

## 参考资料

- [pybinding 官方文档](https://pybinding.site/)
- Kekulé 与高阶拓扑绝缘体相关文献

---

*项目已完成。如有问题欢迎提交 Issue。*
