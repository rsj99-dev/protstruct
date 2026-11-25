# protstruct PSIPRED 版

`protstruct` 是一套轻量的蛋白质二级结构与一致性分析工具。本版本采用 **PSIPRED V4.1 单序列模式**，将二级结构预测、二硫键配对、盐桥候选识别、多序列比对、二维一致性评分和二级结构数据库检索整合在同一流程中。

- 版本：`0.3.0`
- Python：3.8+
- Python 核心依赖：无
- 许可证：GPLv3

## 主要功能

- **PSIPRED V4.1 单序列预测**
  - Rust 实现 `seq2mtx` 伪 PSSM 生成。
  - `psipred` 两级神经网络与 `psipass2` 平滑。
  - 输出 H/E/C 三态，β-转角并入 C。
- **FASTA 输入**
  - 支持 `.fa`、`.faa` 和 `.fasta`。
  - 每次比对支持 1–20 条序列。
  - 自动识别 UTF-8、UTF-16 和 GBK 编码。
- **结构特征分析**
  - 动态规划二硫键配对。
  - 基于序列窗口的盐桥候选评分。
  - 可选用 PDB 同链接触进行结构级核验。
- **多序列比对**
  - center-star 策略。
  - Needleman–Wunsch 全局比对，BLOSUM62 打分，仿射空位罚分。
- **一致性评分**
  - 二级结构相同 +1，不同 0。
  - 二硫键/盐桥特征列对匹配 +0.5。
  - 允许跨类型匹配，同一比对列对只计一次。
  - CSV 同时输出总分和总分一致率。
- **数据库检索**
  - PSIPRED H/E/C 预测、k-mer 预筛和 Smith–Waterman 检索。
  - 可选 Rust 并行核心，未构建时回退 Python。
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

如需并行数据库检索，请构建检索核心：

```powershell
cargo build --release --manifest-path rust/protstruct_search/Cargo.toml
```

可通过以下环境变量指定程序路径：

```text
PROTSTRUCT_PSIPRED_BIN
PROTSTRUCT_PLOT_BIN
PROTSTRUCT_SEARCH_BIN
```

## 使用

基本分析：

```powershell
protstruct examples/demo.fasta
```

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

## 测试

```powershell
python -m unittest discover -s tests -v
```

共 37 个测试用例。图相关用例需要 Rust `protstruct-plot` 程序，否则自动跳过；检索核心不可用时，Rust 检索相关用例也会跳过。

## 目录结构

```text
software/
├── protstruct/                  # Python 分析包
├── rust/
│   ├── protstruct_psipred/      # Rust PSIPRED 推理
│   ├── protstruct_plot/         # 零依赖 SVG 绘图
│   └── protstruct_search/       # 可选并行检索核心
├── third_party/psipred/         # PSIPRED V4.1 权重与许可证
├── examples/demo.fasta
├── tests/
├── README.md                    # 英文文档
├── README_zh.md                 # 中文文档
└── CHANGELOG_v0.3.0_vs_v0.2.0.txt
```

## 许可证

PSIPRED 由 David T. Jones 等开发，采用 GPLv3。许可证见 `third_party/psipred/LICENSE`。本软件包整体按 GPLv3 分发。
