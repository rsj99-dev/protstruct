# Model data location

The PSIPRED edition does not load prediction weights from `protstruct/data/`.
The bundled PSIPRED V4.1 weight files are located at:

```text
third_party/psipred/data/
```

The Rust inference core reads these files through the path configured by
`protstruct.psipred_engine`. No additional model download is required.
