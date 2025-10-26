# Kekulé Lattice Simulation

本仓库用于在 **pybinding** 框架下构建 Kekulé-type 蜂窝晶格并计算其能带、波函数与光电导等物理量。代码最初用于探索拓扑角态、边界态对电导的贡献，可用于复现论文中的数值结果或作为进一步研究的起点。

## 功能概览

| 模块 | 作用 | 说明 |
| ---- | ---- | ---- |
| `basic_function.py` | 定义 Kekulé 晶格、裁剪区域、能带/波函数绘制等基础工具。 | 提供 `O_keku`、`trapezoid`、`probability` 等核心函数。
| `Calculation_total.py` | 计算整体光电导。 | 构建模型、求解本征态并对所有能带积分得到 `σ_{xx}`，结果输出到 `result/anti-corner...` 目录。
| `Calculation_each.py` | 分解不同态对光电导的贡献。 | 以角态、边缘态、体态分类累积结果，便于对比分析。
| `EnergyBand_wavefunction.py` | 计算能带并导出波函数概率分布。 | 用于生成 `.eps` 能带图以及波函数 PNG。
| `lattice_graph.py` | 绘制晶格几何示意图。 | 显示子晶格位置与跃迁。
| `gamma.py` / `Mu.py` / `tempreture.py` | 参数扫描脚本。 | 分别扫描展宽 `Γ`、化学势 `μ`、温度 `k_B T` 对结果的影响。
| `opticalConductivity1.py`、`all conductivity.py` | 光电导计算。 | 可作不同跃迁成分之间的比较。
| `wavefuction_test.py`、`compare.py` | 调试脚本。 | 用于测试波函数与其它分析。

## 环境依赖

- Python 3.8.5 （考虑到`pybinding`的不支持过于新的版本）
- [pybinding](https://pybinding.site/) 及其依赖（需要 C++ 编译工具链）
- `numpy`
- `matplotlib`
- `tqdm`
- `numba`

建议使用虚拟环境管理依赖：

```bash
python -m venv .venv
source .venv/bin/activate  # Windows 使用 .venv\\Scripts\\activate
pip install numpy matplotlib tqdm numba
# 安装 pybinding 时请参考官方文档，根据操作系统准备额外依赖
```

## 快速开始

1. **选择脚本**：根据需求运行 `Calculation_total.py`、`Calculation_each.py` 或其它脚本。
2. **运行模拟**：
   ```bash
   python Calculation_total.py
   ```
   脚本会构建六边形裁剪的 Kekulé 晶格、求解本征值/本征矢，并输出能带、模型图以及速度算符矩阵。
3. **查看结果**：
   - 能带与模型图保存在 `result/anti-corner...` 目录下的 `.eps` 文件。
   - 光电导数据写入同目录的 `total_sigma_xx.txt` 等文本文件，可使用 `numpy.loadtxt` 或绘图工具进行分析。

## 配置说明

- 主要参数（如晶格边长 `edgelength`、跃迁振幅 `t1/t2`、温度 `k_B T`、化学势 `Mu`、展宽 `Gamma` 等）在对应脚本顶部集中定义，可根据研究需要修改。
- 若需要批量计算不同参数，可直接调用 `gamma.py`、`Mu.py` 或 `tempreture.py`，这些脚本会循环调用主计算程序并整理结果。
- 输出目录默认为 `result/...`，运行前可手动清理旧数据以避免混淆。

## 可视化与后处理

- 生成的 `.eps` 图像可用 `matplotlib`、`Inkscape` 或其它矢量图工具打开。
- `matrix_output*.txt`、`Vx.txt`、`Vy.txt` 等文本包含速度算符或哈密顿量矩阵，可在 `NumPy`/`Matlab` 中加载继续分析。

## 参考资料

- [pybinding 官方文档](https://pybinding.site/)
- 关于 Kekulé 变形与高阶拓扑相的综述文章

如果在使用过程中遇到问题，欢迎提交 Issue 或 Pull Request 进行讨论与改进。
