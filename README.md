<p align="center">
<<<<<<< Updated upstream
  <img src="docs/logo.png" alt="protstruct logo" width="460">
</p>

<p align="center">
  <img alt="version" src="https://img.shields.io/badge/version-v0.3.0-2563eb">
  <img alt="python" src="https://img.shields.io/badge/python-3.8%2B-3776ab">
  <img alt="license" src="https://img.shields.io/badge/license-MIT-22c55e">
  <a href="https://github.com/rsj99-dev/protstruct"><img alt="GitHub" src="https://img.shields.io/badge/GitHub-rsj99--dev%2Fprotstruct-181717?logo=github"></a>
</p>

# protstruct
=======
  <img src="docs/logo.svg" alt="protstruct logo" width="460">
</p>

# protstruct PSIPRED edition
>>>>>>> Stashed changes

`protstruct` is a lightweight protein secondary-structure and consistency-analysis toolkit. This edition uses **PSIPRED V4.1 in single-sequence mode** and combines secondary-structure prediction, disulfide pairing, salt-bridge candidates, multiple-sequence alignment, two-dimensional consistency scoring, and secondary-structure database search.

- Version: `0.4.0`
- Python: 3.8+
- Core Python dependencies: none
- License: GPLv3

## Features

- **PSIPRED V4.1 single-sequence prediction**
  - Rust implementation of `seq2mtx` pseudo-PSSM generation.
  - Two-stage `psipred` neural network and `psipass2` smoothing.
  - Batch FASTA inference loads the weights once for all sequences.
  - H/E/C output; β-turns are included in C.
- **FASTA input**
  - `.fa`, `.faa`, and `.fasta`.
  - No sequence-count limit per comparison.
  - UTF-8, UTF-16, and GBK decoding.
- **Structural features**
  - Batch Rust dynamic-programming disulfide pairing with a Python fallback.
  - Sequence-window salt-bridge candidate scoring.
  - Optional same-chain PDB contact validation.
- **Multiple-sequence alignment**
  - Center-star strategy.
  - Batch Rust Needleman–Wunsch with BLOSUM62 and affine gaps.
  - Configurable 1–64 threads; `0` selects all visible logical CPUs, capped at 64.
  - Matching pure-Python fallback when the Rust core is unavailable.
- **Consistency scoring**
  - Secondary-structure match: +1; mismatch: 0.
  - Disulfide/salt-bridge feature-pair match: +0.5.
  - Cross-type feature matches are allowed and counted once per aligned pair.
  - CSV output includes both total score and score ratio.
- **Database search**
  - Full Rust pipeline: batch PSIPRED prediction, k-mer prefilter, parallel
    Smith–Waterman search, filtering, and ranking.
  - Pure-Python fallback with equivalent results.
  - When a cache path is supplied, Python manages cached predictions and
    invokes the Rust candidate scorer.
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

To enable the full Rust database-search pipeline, build the search core:

```powershell
cargo build --release --manifest-path rust/protstruct_search/Cargo.toml
```

To accelerate pairwise alignment and disulfide pairing, build the batch core:

```powershell
cargo build --release --manifest-path rust/protstruct_core/Cargo.toml
```

Custom binary paths can be set with:

```text
PROTSTRUCT_PSIPRED_BIN
PROTSTRUCT_PLOT_BIN
PROTSTRUCT_SEARCH_BIN
PROTSTRUCT_CORE_BIN
```

## Usage

Basic analysis:

```powershell
protstruct examples/demo.fasta
```

Use an explicit thread count:

```powershell
protstruct examples/demo.fasta --threads 16
```

The standalone GUI exposes the same setting in its **Threads** field. The value
is clamped to 1–64; the active count is further limited by visible logical CPUs
and available alignment tasks. Rust alignment parallelizes center-to-sequence
pairs and disulfide pairing. The GUI also predicts analysis sequences in
parallel and parallelizes database-search prediction and Smith-Waterman scoring
through the dedicated **Threads** field on the search page.

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

## Threading benchmark

```powershell
python tests/benchmark_threads.py --sequences 16 --length 1500
```

This compares the original pure-Python path, Rust with 1/8/16/64 requested
threads, and automatic thread selection. Peak memory includes the Python parent
and Rust core subprocess.

## Tests

```powershell
python -m unittest discover -s tests -v
```

The suite contains 45 tests. Figure tests require the Rust `protstruct-plot` binary; otherwise they are skipped. Rust core/search tests are optional when their binaries are unavailable.

## Directory layout

```text
software/
├── protstruct/                  # Python analysis package
├── rust/
│   ├── protstruct_psipred/      # Rust PSIPRED inference
│   ├── protstruct_plot/         # Dependency-free SVG plotting
│   ├── protstruct_search/       # Full Rust prediction and search pipeline
│   └── protstruct_core/         # Batch alignment and disulfide core
├── third_party/psipred/         # PSIPRED V4.1 weights and license
├── examples/demo.fasta
├── tests/
├── README.md                    # English documentation
├── README_zh.md                 # Chinese documentation
├── CHANGELOG_v0.3.0_vs_v0.2.0.txt
└── CHANGELOG_v0.4.0_vs_v0.3.0.txt
```

## License

PSIPRED was developed by David T. Jones and colleagues and is distributed under GPLv3. The PSIPRED license is included at `third_party/psipred/LICENSE`. This software package is distributed under GPLv3.
