"""Rust SVG 绘图入口（无 matplotlib / numpy 依赖）。

Python 只负责把分析结果写成简单的 TSV，然后调用
``software/rust/protstruct_plot`` 编译出的 ``protstruct-plot`` 生成 SVG。
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from .analysis import AnalysisResult
from .consistency import ConsistencyReport, compute_consistency


def _plot_binary() -> str:
    configured = os.environ.get("PROTSTRUCT_PLOT_BIN")
    if configured and Path(configured).is_file():
        return configured
    root = Path(__file__).resolve().parent.parent
    candidates = [
        root / "rust" / "protstruct_plot" / "target" / "release" / "protstruct-plot",
        root / "rust" / "protstruct_plot" / "target" / "release" / "protstruct-plot.exe",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    raise RuntimeError(
        "Rust SVG plotter protstruct-plot was not found. Build it first:\n"
        "  cd software/rust/protstruct_plot && cargo build --release\n"
        "or set PROTSTRUCT_PLOT_BIN."
    )


def _run_plot(mode: str, output: str, data: Optional[str] = None, lang: str = "en") -> str:
    path = Path(output)
    if path.suffix.lower() != ".svg":
        raise ValueError("The Rust SVG backend only supports .svg output.")
    path.parent.mkdir(parents=True, exist_ok=True)
    input_path = None
    try:
        command = [_plot_binary(), mode, "--output", str(path), "--lang", lang]
        if data is not None:
            with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False, encoding="utf-8") as handle:
                handle.write(data)
                input_path = handle.name
            command.extend(["--input", input_path])
        process = subprocess.run(command, check=True, capture_output=True, text=True, encoding="utf-8")
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError("The Rust SVG plotter did not create a file: %s" % path)
        return str(path)
    except FileNotFoundError as exc:
        raise RuntimeError("protstruct-plot was not found; build the Rust plotting core first.") from exc
    except subprocess.CalledProcessError as exc:
        message = (exc.stderr or exc.stdout or "").strip()
        raise RuntimeError("protstruct-plot failed: %s" % message) from exc
    finally:
        if input_path:
            try:
                os.remove(input_path)
            except OSError:
                pass


def plot_workflow(path: str, lang: str = "en") -> str:
    """Generate the workflow SVG."""
    return _run_plot("workflow", path, lang=lang)


def _consistency_payload(result: AnalysisResult, report: ConsistencyReport) -> str:
    lines = [
        "REF\t%s" % report.reference,
        "MAX\t%.6f" % report.max_total,
        "WIDTH\t%d" % report.width,
        "AGREE\t%s" % ",".join("%.6f" % value for value in report.agreement),
    ]
    for score, ss in zip(report.scores, report.aligned_ss):
        lines.append("SS\t%s\t%s" % (score.identifier, ss))
    for score in report.scores:
        lines.append("SCORE\t%s\t%.6f\t%d" % (score.identifier, score.total, 1 if score.is_reference else 0))
    for row, score in enumerate(report.scores):
        for columns in score.disulfide_columns:
            for column in columns:
                lines.append("DS\t%d\t%d" % (row, column))
        for columns in score.salt_bridge_columns:
            for column in columns:
                lines.append("SB\t%d\t%d" % (row, column))
    return "\n".join(lines) + "\n"


def plot_consistency(
    result: AnalysisResult,
    report: ConsistencyReport,
    path: str,
    dpi: int = 200,
    figure_scale: float = 1.0,
    lang: str = "en",
) -> str:
    """Generate the consistency SVG.

    dpi/figure_scale are retained for compatibility; SVG is vector output and
    does not require DPI.
    """
    del dpi, figure_scale
    return _run_plot("consistency", path, _consistency_payload(result, report), lang=lang)


def _short_label(identifier: str) -> str:
    return identifier.split("（")[0].split("(")[0].strip()


def _mtba_payload(
    paralog_result: AnalysisResult,
    ortholog_result: AnalysisResult,
    paralog_labels: Optional[Sequence[str]],
    ortholog_labels: Optional[Sequence[str]],
) -> str:
    paralog_report = compute_consistency(paralog_result)
    ortholog_report = compute_consistency(ortholog_result)
    labels_paralog = list(paralog_labels) if paralog_labels else [
        _short_label(item.identifier) for item in paralog_result.items
    ]
    labels_ortholog = list(ortholog_labels) if ortholog_labels else [
        item.identifier for item in ortholog_result.items
    ]
    lines = []
    for label, item, score in zip(labels_paralog, paralog_result.items, paralog_report.scores):
        fractions = item.fractions()
        lines.append(
            "PARA\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
            % (label, fractions["H"] * 100, fractions["E"] * 100, fractions["C"] * 100, score.total, score.score_rate)
        )
    for label, item, score in zip(labels_ortholog, ortholog_result.items, ortholog_report.scores):
        fractions = item.fractions()
        lines.append(
            "ORTH\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
            % (label, fractions["H"] * 100, fractions["E"] * 100, fractions["C"] * 100, score.total, score.score_rate)
        )
    return "\n".join(lines) + "\n"


def plot_mtba(
    paralog_result: AnalysisResult,
    ortholog_result: AnalysisResult,
    path: str,
    paralog_labels: Optional[Sequence[str]] = None,
    ortholog_labels: Optional[Sequence[str]] = None,
    lang: str = "en",
) -> str:
    """Generate the combined composition-and-score SVG."""
    return _run_plot(
        "mtba",
        path,
        _mtba_payload(paralog_result, ortholog_result, paralog_labels, ortholog_labels),
        lang=lang,
    )
