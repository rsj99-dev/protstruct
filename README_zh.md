<p align="center">
  <img src="docs/logo.svg" alt="protstruct logo" width="460">
</p>

# protstruct PSIPRED 版

`protstruct` 是一套轻量的蛋白质二级结构与一致性分析工具。本版本采用 **PSIPRED V4.1 单序列模式**，将二级结构预测、二硫键配对、盐桥候选识别、多序列比对、二维一致性评分和二级结构数据库检索整合在同一流程中。

- 版本：`0.4.0`
- Python：3.8+
- Python 核心依赖：无
- 许可证：GPLv3

## 主要功能

- **PSIPRED V4.1 单序列预测**
  - Rust 实现 `seq2mtx` 伪 PSSM 生成。
  - `psipred` 两级神经网络与 `psipass2` 平滑。
  - 支持批量 FASTA 推理，整批序列只加载一次权重。
  - 输出 H/E/C 三态，β-转角并入 C。
- **FASTA 输入**
  - 支持 `.fa`、`.faa` 和 `.fasta`。
  - 每次比对的序列数量不限。
  - 自动识别 UTF-8、UTF-16 和 GBK 编码。
- **结构特征分析**
  - Rust 批量动态规划二硫键配对，保留 Python 回退。
  - 基于序列窗口的盐桥候选评分。
  - 可选用 PDB 同链接触进行结构级核验。
- **多序列比对**
  - center-star 策略。
  - Rust 批量 Needleman–Wunsch 全局比对，BLOSUM62 打分，仿射空位罚分。
  - 支持 1–64 线程；`0` 表示使用全部可见逻辑核心，最多 64 个。
  - Rust 核心不可用时自动使用结果一致的 Python 回退。
- **一致性评分**
  - 二级结构相同 +1，不同 0。
  - 二硫键/盐桥特征列对匹配 +0.5。
  - 允许跨类型匹配，同一比对列对只计一次。
  - CSV 同时输出总分和总分一致率。
- **数据库检索**
  - PSIPRED H/E/C 预测、k-mer 预筛和 Smith–Waterman 检索。
  - Rust 核心完成批量预测、k-mer 预筛、并行比对、过滤和排序。
  - 未构建检索核心时回退 Python。
  - 使用缓存路径时，由 Python 管理预测缓存并调用 Rust 候选评分。
- **英文输出**
  - CLI、报告、CSV、检索 TSV 和默认 SVG 标签均为英文。
- **仅输出 SVG**
  - Rust 绘图后端只生成 SVG。
  - 不依赖 matplotlib 或 numpy。

## 安装

在当前目录执行：

```powershell
pip install -e .
```

PSIPRED V4.1 权重位于：

```text
third_party/psipred/data/
```

### Rust 组件

Windows 预编译推理程序位于：

```text
rust/protstruct_psipred/bin/win/protstruct-psipred.exe
```

Linux 或 macOS 请构建 PSIPRED 推理核心：

```powershell
cargo build --release --manifest-path rust/protstruct_psipred/Cargo.toml
```

如需生成 SVG，请构建绘图核心：

```powershell
cargo build --release --manifest-path rust/protstruct_plot/Cargo.toml
```

如需完整的 Rust 数据库检索流程，请构建检索核心：

```powershell
cargo build --release --manifest-path rust/protstruct_search/Cargo.toml
```

如需加速两两比对和二硫键配对，请构建批量核心：

```powershell
cargo build --release --manifest-path rust/protstruct_core/Cargo.toml
```

可通过以下环境变量指定程序路径：

```text
PROTSTRUCT_PSIPRED_BIN
PROTSTRUCT_PLOT_BIN
PROTSTRUCT_SEARCH_BIN
PROTSTRUCT_CORE_BIN
```

## 使用

基本分析：

```powershell
protstruct examples/demo.fasta
```

指定线程数：

```powershell
protstruct examples/demo.fasta --threads 16
```

独立 GUI 的 **Threads** 输入框提供相同设置，数值限制在 1–64；实际线程数还会受可见逻辑核心数和可用比对任务数限制。Rust 比对会在中心序列与其余序列之间并行，并批量并行二硫键配对；GUI 还会并行执行分析序列的 PSIPRED 预测，并通过检索页独立的 **Threads** 输入框并行执行数据库预测与 Smith–Waterman 评分。

导出报告：

```powershell
protstruct examples/demo.fasta -f json -o report.json
```

生成一致性 SVG 和评分 CSV：

```powershell
protstruct examples/demo.fasta --figure consistency.svg --scores scores.csv
```

检索序列数据库：

```powershell
protstruct-search query.fasta database.fasta --top 20 -o hits.tsv
```

`--figure` 只接受 `.svg`。只提供 `--figure` 或 `--scores` 时，程序会自动推导另一个输出路径。

## Python API

```python
from protstruct import analyze_fasta, compute_consistency
from protstruct.figure import plot_consistency

result = analyze_fasta("examples/demo.fasta")
report = compute_consistency(result)
plot_consistency(result, report, "consistency.svg")
```

## PSIPRED 使用范围

本包实现 PSIPRED V4.1 单序列模式，不运行 PSI-BLAST，也不构建 UniRef90 profile，因此预测结果不等同于完整 profile 模式。包内不包含 C/C++ 实现，只在 `third_party/psipred/` 中保留 PSIPRED V4.1 权重。

短肽、低复杂度序列、膜蛋白和富含二硫键的蛋白可能出现较大偏差。预测的二硫键和盐桥应结合三维结构或实验证据判断。

## 线程基准

```powershell
python tests/benchmark_threads.py --sequences 16 --length 1500
```

该脚本比较原始 Python、Rust 1/8/16/64 请求线程和自动线程模式，内存峰值包含
Python 父进程与 Rust 核心子进程。

## 测试

```powershell
python -m unittest discover -s tests -v
```

共 45 个测试用例。图相关用例需要 Rust `protstruct-plot` 程序，否则自动跳过；批量核心或检索核心不可用时，对应 Rust 用例也会跳过。

## 目录结构

```text
software/
├── protstruct/                  # Python 分析包
├── rust/
│   ├── protstruct_psipred/      # Rust PSIPRED 推理
│   ├── protstruct_plot/         # 零依赖 SVG 绘图
│   ├── protstruct_search/       # Rust 批量预测与完整检索
│   └── protstruct_core/         # Rust 批量比对与二硫键核心
├── third_party/psipred/         # PSIPRED V4.1 权重与许可证
├── examples/demo.fasta
├── tests/
├── README.md                    # 英文文档
├── README_zh.md                 # 中文文档
├── CHANGELOG_v0.3.0_vs_v0.2.0.txt
└── CHANGELOG_v0.4.0_vs_v0.3.0.txt
```

## 许可证

PSIPRED 由 David T. Jones 等开发，采用 GPLv3。许可证见 `third_party/psipred/LICENSE`。本软件包整体按 GPLv3 分发。
