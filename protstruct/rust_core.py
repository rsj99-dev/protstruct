"""Rust 批量计算核心的通用调用入口。"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

MAX_THREADS = 64


def effective_threads(threads: Optional[int] = 0) -> int:
    """返回 1–64 的实际线程数；0 表示使用可见逻辑核心（最多 64）。"""
    available = os.cpu_count() or 1
    requested = available if not threads else threads
    return max(1, min(int(requested), available, MAX_THREADS))


def _core_binary() -> Optional[str]:
    """查找批量核心，支持环境变量和本地 release 构建路径。"""
    configured = os.environ.get("PROTSTRUCT_CORE_BIN")
    if configured and Path(configured).is_file():
        return configured
    root = Path(__file__).resolve().parent.parent
    for name in ("protstruct-core.exe", "protstruct-core"):
        for directory in (
            root / "rust" / "protstruct_core" / "target" / "release",
            root / "rust" / "protstruct_core" / "bin" / "win",
        ):
            candidate = directory / name
            if candidate.is_file():
                return str(candidate)
    return None


def _run_core(mode: str, payload: str, threads: Optional[int] = 0) -> Optional[str]:
    """运行批量核心；任何缺失、协议或执行错误都交由 Python 回退。"""
    binary = _core_binary()
    if not binary:
        return None
    try:
        process = subprocess.run(
            [binary, mode, "--threads", str(effective_threads(threads))],
            input=payload,
            capture_output=True,
            check=True,
            text=True,
            encoding="utf-8",
        )
        return process.stdout
    except (OSError, subprocess.SubprocessError, UnicodeError):
        return None


def align_pairs_rust(
    pairs: Sequence[Tuple[str, str, str]],
    gap_open: float = -11.0,
    gap_extend: float = -1.0,
    threads: Optional[int] = 0,
) -> Optional[List[Tuple[str, str, float]]]:
    """批量两两比对，返回顺序与输入一致的比对字符串和得分。"""
    if not pairs:
        return []
    lines = ["%.17g\t%.17g\t0" % (gap_open, gap_extend)]
    lines.extend("%d\t%s\t%s" % (index, first, second) for index, (_, first, second) in enumerate(pairs))
    output = _run_core("align", "\n".join(lines) + "\n", threads)
    if output is None:
        return None
    try:
        parsed = {}
        for line in output.splitlines():
            fields = line.split("\t")
            if len(fields) != 4:
                return None
            index = int(fields[0])
            if index in parsed:
                return None
            parsed[index] = (fields[1], fields[2], float(fields[3]))
        if set(parsed) != set(range(len(pairs))):
            return None
        results = []
        for index, (_, first, second) in enumerate(pairs):
            left, right, score = parsed[index]
            if left.replace("-", "") != first or right.replace("-", "") != second or len(left) != len(right):
                return None
            results.append((left, right, score))
        return results
    except (ValueError, IndexError):
        return None


def disulfide_pairs_rust(
    records: Sequence[Tuple[str, str, str, int, int]],
    threads: Optional[int] = 0,
) -> Optional[List[List[Tuple[int, int]]]]:
    """批量二硫键配对，记录格式为 sequence、ss、最小间隔、最大环长。"""
    if not records:
        return []
    lines = [
        "%d\t%s\t%s\t%d\t%d" % (index, sequence, ss or "-", min_separation, max_loop)
        for index, (sequence, ss, min_separation, max_loop) in enumerate(records)
    ]
    output = _run_core("disulfide", "\n".join(lines) + "\n", threads)
    if output is None:
        return None
    pairs: List[List[Tuple[int, int]]] = [[] for _ in records]
    try:
        for line in output.splitlines():
            fields = line.split("\t")
            index = int(fields[0])
            if index < 0 or index >= len(records) or len(fields) not in (1, 3):
                return None
            if len(fields) == 3:
                first, second = int(fields[1]), int(fields[2])
                sequence, _, min_separation, _ = records[index]
                sequence = sequence.upper()
                if not (1 <= first < second <= len(sequence)):
                    return None
                if sequence[first - 1] != "C" or sequence[second - 1] != "C" or second - first < min_separation:
                    return None
                pairs[index].append((first, second))
        for record_pairs in pairs:
            used = set()
            for first, second in record_pairs:
                if first in used or second in used:
                    return None
                used.update((first, second))
        return pairs
    except (ValueError, IndexError):
        return None
