"""多序列比对线程扩展性与峰值内存基准。

用法：
    python tests/benchmark_threads.py
    python tests/benchmark_threads.py --sequences 32 --length 700
"""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

try:
    import psutil
except ImportError:  # 仅影响内存统计
    psutil = None


def make_sequences(count: int, length: int, identity: float = 0.9) -> List[str]:
    """生成同一家族的突变序列，避免随机序列造成异常宽的多序列比对。"""
    random.seed(20261125)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    base = [random.choice(alphabet) for _ in range(length)]
    sequences = []
    for _ in range(count):
        sequence = list(base)
        for position in range(length):
            if random.random() >= identity:
                sequence[position] = random.choice(alphabet)
        sequences.append("".join(sequence))
    return sequences


def run_worker(mode: str, sequences: int, length: int, identity: float) -> Dict[str, object]:
    from unittest import mock

    from protstruct.alignment import build_alignment
    from protstruct.rust_core import effective_threads

    data = make_sequences(sequences, length, identity)
    identifiers = ["seq%d" % index for index in range(sequences)]
    started = time.perf_counter()
    if mode == "python":
        with mock.patch("protstruct.alignment.align_pairs_rust", return_value=None):
            result = build_alignment(identifiers, data, threads=1)
        requested_threads = 1
    else:
        requested_threads = 0 if mode == "auto" else int(mode)
        result = build_alignment(identifiers, data, threads=requested_threads)
    elapsed = time.perf_counter() - started
    return {
        "mode": mode,
        "elapsed": elapsed,
        "requested_threads": requested_threads,
        "threads": 1 if mode == "python" else effective_threads(requested_threads),
        "rows": len(result.aligned),
        "width": result.width,
        "checksum": round(sum(result.identity[0]), 6),
    }


def _tree_rss(process) -> int:
    if psutil is None:
        return 0
    try:
        total = process.memory_info().rss
        for child in process.children(recursive=True):
            try:
                total += child.memory_info().rss
            except psutil.Error:
                pass
        return total
    except psutil.Error:
        return 0


def run_parent(args: argparse.Namespace) -> None:
    # 先测 Rust，避免长时间 Python 基准造成 CPU 降频后影响并行扩展性结果。
    modes = ["1", "8", "16", "64", "auto", "python"]
    rows = []
    for mode in modes:
        command = [
            sys.executable, str(Path(__file__).resolve()),
            "--worker", mode,
            "--sequences", str(args.sequences),
            "--length", str(args.length),
            "--identity", str(args.identity),
        ]
        started = time.perf_counter()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8")
        peak = 0
        monitored = None
        if psutil is not None:
            monitored = psutil.Process(process.pid)
        while process.poll() is None:
            peak = max(peak, _tree_rss(monitored))
            time.sleep(0.005)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise SystemExit("worker %s failed:\n%s\n%s" % (mode, stdout, stderr))
        result = json.loads(stdout.strip().splitlines()[-1])
        result["wall"] = time.perf_counter() - started
        result["peak_mb"] = peak / (1024 * 1024) if peak else None
        rows.append(result)
        print(
            "%-6s requested=%-2s actual=%-2s time=%7.3fs peak=%s rows=%d width=%d"
            % (
                mode,
                result["requested_threads"],
                result["threads"],
                result["elapsed"],
                "%.1f MiB" % result["peak_mb"] if result["peak_mb"] else "n/a",
                result["rows"],
                result["width"],
            )
        )

    python_row = next(row for row in rows if row["mode"] == "python")
    baseline = python_row["elapsed"]
    single = next(row for row in rows if row["mode"] == "1")
    default = next(row for row in rows if row["mode"] == "auto")
    print("\nSummary")
    print("Python original: %.3fs, %.1f MiB" % (baseline, python_row["peak_mb"] or 0.0))
    print("Rust 1 thread  : %.3fs, %.1f MiB, %.1fx faster than Python"
          % (single["elapsed"], single["peak_mb"] or 0.0, baseline / single["elapsed"]))
    print("Rust default   : %d threads actual, %.3fs, %.1f MiB, %.1fx faster than Python, %.1fx faster than Rust-1"
          % (default["threads"], default["elapsed"], default["peak_mb"] or 0.0,
             baseline / default["elapsed"], single["elapsed"] / default["elapsed"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", choices=["python", "auto", "1", "8", "16", "64"])
    parser.add_argument("--sequences", type=int, default=24)
    parser.add_argument("--length", type=int, default=600)
    parser.add_argument("--identity", type=float, default=0.9)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(run_worker(args.worker, args.sequences, args.length, args.identity)))
    else:
        run_parent(args)


if __name__ == "__main__":
    main()
