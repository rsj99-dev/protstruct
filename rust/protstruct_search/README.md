# protstruct Rust search core

Full Rust secondary-structure search pipeline used by:

```powershell
protstruct-search query.fa database.fa --top 20
# Direct Rust invocation
rust/protstruct_search/target/release/protstruct-search --query query.fa --database database.fa --data third_party/psipred/data --top 20
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

In full mode, Rust reads the query and database FASTA files, performs batch
PSIPRED prediction with weights loaded once, applies sequence/SS k-mer
prefiltering, and runs parallel Smith-Waterman scoring with Rayon. Filtering
and final ranking are also performed in Rust.

The legacy `--input` mode still accepts a pre-scored TSV from Python. If the
binary is unavailable, Python uses its equivalent pure-Python pipeline.
