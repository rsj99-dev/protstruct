# protstruct-plot

Dependency-free Rust SVG generator for the paper figures.

```powershell
cd software/rust/protstruct_plot
cargo build --release
```

Modes:

```text
protstruct-plot workflow --output fig1.svg
protstruct-plot consistency --input consistency.tsv --output fig2.svg
protstruct-plot mtba --input mtba.tsv --output fig3.svg
```

No Python plotting dependency is used; the SVG text references
`Times New Roman / SimSun` font fallbacks. Python calls this binary from
`protstruct.figure` and the paper reproduction scripts.
