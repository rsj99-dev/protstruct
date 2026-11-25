<p align="center">
  <img src="docs/logo.png" alt="protstruct logo" width="460">
</p>

<p align="center">
  <img alt="version" src="https://img.shields.io/badge/version-v0.3.0-2563eb">
  <img alt="python" src="https://img.shields.io/badge/python-3.8%2B-3776ab">
  <img alt="license" src="https://img.shields.io/badge/license-MIT-22c55e">
  <a href="https://github.com/rsj99-dev/protstruct"><img alt="GitHub" src="https://img.shields.io/badge/GitHub-rsj99--dev%2Fprotstruct-181717?logo=github"></a>
</p>

# protstruct PSIPRED edition

`protstruct` is a lightweight protein secondary-structure and consistency-analysis toolkit. This edition uses **PSIPRED V4.1 in single-sequence mode** and combines secondary-structure prediction, disulfide pairing, salt-bridge candidates, multiple-sequence alignment, two-dimensional consistency scoring, and secondary-structure database search.

- Version: `0.3.0`
- Python: 3.8+
- Core Python dependencies: none
- License: GPLv3

## Features

- **PSIPRED V4.1 single-sequence prediction**
  - Rust implementation of `seq2mtx` pseudo-PSSM generation.
  - Two-stage `psipred` neural network and `psipass2` smoothing.
  - H/E/C output; β-turns are included in C.
- **FASTA input**
  - `.fa`, `.faa`, and `.fasta`.
  - 1–20 sequences per comparison.
  - UTF-8, UTF-16, and GBK decoding.
- **Structural features**
  - Dynamic-programming disulfide pairing.
  - Sequence-window salt-bridge candidate scoring.
  - Optional same-chain PDB contact validation.
- **Multiple-sequence alignment**
  - Center-star strategy.
  - Needleman–Wunsch with BLOSUM62 and affine gaps.
- **Consistency scoring**
  - Secondary-structure match: +1; mismatch: 0.
  - Disulfide/salt-bridge feature-pair match: +0.5.
  - Cross-type feature matches are allowed and counted once per aligned pair.
  - CSV output includes both total score and score ratio.
- **Database search**
  - PSIPRED H/E/C prediction, k-mer prefilter, and Smith–Waterman search.
  - Optional Rust parallel core with a Python fallback.
- **English output**
  - CLI messages, reports, CSV, search TSV, and default SVG labels are in English.
- **SVG-only plotting**
  - The Rust plotting backend emits SVG only.
  - No matplotlib or numpy dependency.

## Installation

From this directory:

```powershell
pip install -e .
```

The PSIPRED V4.1 weights are included under:

```text
third_party/psipred/data/
```

### Rust components

The Windows PSIPRED binary is included at:

```text
rust/protstruct_psipred/bin/win/protstruct-psipred.exe
```

On Linux or macOS, build the PSIPRED core:

```powershell
cargo build --release --manifest-path rust/protstruct_psipred/Cargo.toml
```

To enable SVG plotting, build the plotting core:

```powershell
cargo build --release --manifest-path rust/protstruct_plot/Cargo.toml
```

To enable parallel database search, build the search core:

```powershell
cargo build --release --manifest-path rust/protstruct_search/Cargo.toml
```

Custom binary paths can be set with:

```text
PROTSTRUCT_PSIPRED_BIN
PROTSTRUCT_PLOT_BIN
PROTSTRUCT_SEARCH_BIN
```

## Usage

Basic analysis:

```powershell
protstruct examples/demo.fasta
```

Export a report:

```powershell
protstruct examples/demo.fasta -f json -o report.json
```

Generate a consistency SVG and score CSV:

```powershell
protstruct examples/demo.fasta --figure consistency.svg --scores scores.csv
```

Search a sequence database:

```powershell
protstruct-search query.fasta database.fasta --top 20 -o hits.tsv
```

`--figure` accepts `.svg` only. When only `--figure` or `--scores` is supplied, the other output path is derived automatically.

## Python API

```python
from protstruct import analyze_fasta, compute_consistency
from protstruct.figure import plot_consistency

result = analyze_fasta("examples/demo.fasta")
report = compute_consistency(result)
plot_consistency(result, report, "consistency.svg")
```

## PSIPRED scope

This package implements PSIPRED V4.1 single-sequence mode. It does not run PSI-BLAST or build a UniRef90 profile, so its predictions are not equivalent to the full profile-based PSIPRED server. The package contains no C/C++ implementation; only the PSIPRED V4.1 weights are retained under `third_party/psipred/`.

Single-sequence predictions can be unreliable for short peptides, low-complexity sequences, membrane proteins, and disulfide-rich proteins. Predicted disulfides and salt bridges should be checked against structural or experimental evidence.

## Tests

```powershell
python -m unittest discover -s tests -v
```

The suite contains 37 tests. Figure tests require the Rust `protstruct-plot` binary; otherwise they are skipped. Rust search tests are optional when the search binary is unavailable.

## Directory layout

```text
software/
├── protstruct/                  # Python analysis package
├── rust/
│   ├── protstruct_psipred/      # Rust PSIPRED inference
│   ├── protstruct_plot/         # Dependency-free SVG plotting
│   └── protstruct_search/       # Optional parallel search core
├── third_party/psipred/         # PSIPRED V4.1 weights and license
├── examples/demo.fasta
├── tests/
├── README.md                    # English documentation
├── README_zh.md                 # Chinese documentation
└── CHANGELOG_v0.3.0_vs_v0.2.0.txt
```

## License

PSIPRED was developed by David T. Jones and colleagues and is distributed under GPLv3. The PSIPRED license is included at `third_party/psipred/LICENSE`. This software package is distributed under GPLv3.
