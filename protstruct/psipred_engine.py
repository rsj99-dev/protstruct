"""Rust PSIPRED V4.1 single-sequence inference bridge.

The neural-network inference (seq2mtx pseudo-PSSM + psipred + psipass2)
is implemented in ``software/rust/protstruct_psipred``.  Only the original
PSIPRED V4.1 weight files remain under ``software/third_party/psipred/data``.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from functools import lru_cache
from typing import Iterable, List
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
        root / "rust" / "protstruct_psipred" / "target" / "release" / "protstruct-psipred.exe",
        root / "rust" / "protstruct_psipred" / "target" / "release" / "protstruct-psipred",
        root / "rust" / "protstruct_psipred" / "bin" / "win" / "protstruct-psipred.exe",
        root / "rust" / "protstruct_psipred" / "bin" / "linux" / "protstruct-psipred",
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


def predict_psipred_batch(sequences: Iterable[str]) -> List[str]:
    """一次启动 Rust 并加载权重，批量预测 FASTA 序列。"""
    residues = [sequence.upper() for sequence in sequences]
    if not residues:
        return []
    if len(residues) == 1:
        return [predict_psipred(residues[0])]
    for sequence in residues:
        if len(sequence) < 5:
            raise ValueError("PSIPRED V4.1 requires a sequence of at least 5 residues.")
        if len(sequence) > 10000:
            raise ValueError("Sequence length exceeds the PSIPRED limit of 10000 residues.")

    path = None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False, encoding="utf-8") as handle:
            path = handle.name
            for index, sequence in enumerate(residues):
                handle.write(">seq%d\n%s\n" % (index, sequence))
        process = subprocess.run(
            [str(_rust_binary()), "--data", str(WEIGHTS_DIR), "--input", path],
            capture_output=True,
            check=False,
        )
        if process.returncode != 0:
            message = (process.stderr or process.stdout or b"").decode("utf-8", "replace").strip()
            # 兼容不支持批量输入参数的旧版单序列二进制。
            if "unknown" in message.lower() or "usage:" in message.lower():
                return [_predict_cached(sequence) for sequence in residues]
            raise RuntimeError("Rust PSIPRED batch inference failed: %s" % message)
        states: List[str] = []
        for line in process.stdout.decode("utf-8", "replace").splitlines():
            value = line.strip()
            if value and not value.startswith(">"):
                states.append(value)
        if len(states) != len(residues):
            # 旧版二进制不支持 --input，可能返回单条状态而非报错。
            return [_predict_cached(sequence) for sequence in residues]
        for sequence, state in zip(residues, states):
            if len(state) != len(sequence) or not set(state) <= {"H", "E", "C"}:
                raise RuntimeError("Rust PSIPRED returned an invalid state string.")
        return states
    finally:
        if path:
            try:
                os.remove(path)
            except OSError:
                pass
