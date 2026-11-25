"""Rust PSIPRED V4.1 single-sequence inference bridge.

The neural-network inference (seq2mtx pseudo-PSSM + psipred + psipass2)
is implemented in ``software/rust/protstruct_psipred``.  Only the original
PSIPRED V4.1 weight files remain under ``software/third_party/psipred/data``.
"""

from __future__ import annotations

import os
import subprocess
from functools import lru_cache
from pathlib import Path

PSIPRED_ROOT = Path(__file__).resolve().parent.parent / "third_party" / "psipred"
WEIGHTS_DIR = PSIPRED_ROOT / "data"


def _rust_binary() -> Path:
    configured = os.environ.get("PROTSTRUCT_PSIPRED_BIN")
    if configured:
        candidate = Path(configured)
        if candidate.is_file():
            return candidate
    root = Path(__file__).resolve().parent.parent
    candidates = [
        root / "rust" / "protstruct_psipred" / "bin" / "win" / "protstruct-psipred.exe",
        root / "rust" / "protstruct_psipred" / "bin" / "linux" / "protstruct-psipred",
        root / "rust" / "protstruct_psipred" / "target" / "release" / "protstruct-psipred.exe",
        root / "rust" / "protstruct_psipred" / "target" / "release" / "protstruct-psipred",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise RuntimeError(
        "Rust PSIPRED binary protstruct-psipred was not found. Build it first:\n"
        "  cd software/rust/protstruct_psipred && cargo build --release\n"
        "or set PROTSTRUCT_PSIPRED_BIN."
    )


@lru_cache(maxsize=512)
def _predict_cached(sequence: str) -> str:
    if len(sequence) < 5:
        raise ValueError("PSIPRED V4.1 requires a sequence of at least 5 residues.")
    if len(sequence) > 10000:
        raise ValueError("Sequence length exceeds the PSIPRED limit of 10000 residues.")
    command = [
        str(_rust_binary()),
        "--data",
        str(WEIGHTS_DIR),
        "--sequence",
        sequence,
    ]
    process = subprocess.run(command, capture_output=True, check=False)
    if process.returncode != 0:
        message = (process.stderr or process.stdout or b"").decode("utf-8", "replace").strip()
        raise RuntimeError("Rust PSIPRED inference failed: %s" % message)
    state = process.stdout.decode("utf-8", "replace").strip()
    if len(state) != len(sequence) or not set(state) <= {"H", "E", "C"}:
        raise RuntimeError("Rust PSIPRED returned an invalid state string.")
    return state


def predict_psipred(sequence: str) -> str:
    """Run the Rust PSIPRED V4.1 single-sequence implementation."""
    return _predict_cached(sequence.upper())
