# protstruct-psipred

Dependency-free Rust implementation of the PSIPRED V4.1 single-sequence
inference path:

```text
seq2mtx pseudo-PSSM
→ three averaged first-level neural networks
→ PSIPASS2 second-level network and filtering
→ H/E/C
```

No C/C++ code or external Rust crates are used. The original PSIPRED V4.1
weight files are in `../../third_party/psipred/data/`.

## Build

```powershell
cd software/rust/protstruct_psipred
cargo build --release
```

The prebuilt Windows binary is in `bin/win/protstruct-psipred.exe`.
Linux/macOS users can build from this source.

## Usage

```powershell
protstruct-psipred --data ../../third_party/psipred/data --sequence ACDEFGHIKL
```

The Python package calls this binary through `protstruct.psipred_engine`.
The implementation is validated against the original PSIPRED V4.1 C
output for the demo and MtbA datasets (identical residue states).

PSIPRED and its weights are GPLv3; see `../../third_party/psipred/LICENSE`.
