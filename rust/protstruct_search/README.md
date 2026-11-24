# protstruct Rust search core

Parallel secondary-structure Smith-Waterman scorer used by:

```powershell
protstruct-search query.fa database.fa --top 20
```

## Build

Requires a Rust toolchain (cargo/rustc).

```powershell
cd software/rust/protstruct_search
cargo build --release
```

On Windows without the MSVC linker (`link.exe`), use:

```powershell
cargo +stable-x86_64-pc-windows-gnu build --release
```

The binary is written to `target/release/protstruct-search` (Windows:
`protstruct-search.exe`). The Python CLI automatically looks for this path;
you can also set:

```powershell
$env:PROTSTRUCT_SEARCH_BIN = "D:\path\to\protstruct-search.exe"
```

If the binary is absent, `protstruct-search` uses the pure-Python fallback.

The Rust core receives a pre-filtered TSV from Python and performs the
secondary-structure local alignment in parallel with Rayon. Sequence k-mer
Jaccard similarity is computed in Rust for tie-breaking.
